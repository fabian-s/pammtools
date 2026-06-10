# ===========================================================================
# IC benchmark: method fitting + estimand extraction
# ---------------------------------------------------------------------------
# One fit_<method>() per method, each returning the standardized long tibble
#   (estimand, t, x, est, se, lower, upper, ambiguity)
# The caller (run-task.R) adds method/cell/rep metadata and joins truth.
#
# CI constructions are deliberately commensurable: all simulation-based CIs
# are empirical quantiles of NSIM coefficient draws pushed through the QOI
# (pamm_ic pooled draws via the package machinery; midpoint/oracle via
# ci_type = "sim"; ic_par via manual normal draws through Weibull closed
# forms). se is reported as the sd of QOI draws where we control the draws,
# else as implied-normal width / (2 * z_{0.975}).
# Requires config.R + dgp.R sourced first; package pammtools loaded.
# ===========================================================================

suppressMessages({
  library(survival)
  library(mgcv)
})

Z975 <- stats::qnorm(0.975)

# vendored from pammtools:::rubin_inflate_qoi_draws (internal; keep in sync)
rubin_inflate_local <- function(draws, estimates) {
  m <- length(draws)
  if (m <= 1L) {
    return(draws)
  }
  qhat <- do.call(cbind, estimates)
  qbar <- rowMeans(qhat)
  between_scale <- sqrt((m + 1) / (m - 1))
  lapply(seq_len(m), function(i) {
    shift <- (between_scale - 1) * (qhat[, i] - qbar)
    sweep(draws[[i]], 1L, shift, `+`)
  })
}

# Rubin's rules with Barnard-Rubin small-sample df for a scalar.
# nu_com: complete-data residual df (large for PED pseudo-data -> close to
# classic Rubin, but reported for transparency).
rubin_scalar_br <- function(est, var, nu_com, level = 0.95) {
  m <- length(est)
  qbar <- mean(est)
  w <- mean(var)
  b <- stats::var(est)
  tot <- w + (1 + 1 / m) * b
  riv <- (1 + 1 / m) * b / w
  df_old <- if (riv <= 0) 1e8 else (m - 1) * (1 + 1 / riv)^2
  gamma <- (1 + 1 / m) * b / tot
  df_obs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - gamma)
  df_br <- 1 / (1 / df_old + 1 / df_obs)
  tc <- stats::qt(1 - (1 - level) / 2, df_br)
  se <- sqrt(tot)
  tibble(
    est = qbar,
    se = se,
    lower = qbar - tc * se,
    upper = qbar + tc * se,
    df = df_br
  )
}

qoi_rows <- function(
  estimand,
  t,
  x,
  est,
  se,
  lower,
  upper,
  ambiguity = NA_real_
) {
  tibble(
    estimand = estimand,
    t = t,
    x = x,
    est = est,
    se = se,
    lower = lower,
    upper = upper,
    ambiguity = ambiguity
  )
}

# ---- prediction grids --------------------------------------------------------
# Full CUT grid per covariate stratum, grouped (survival cumulates within group)
grid_binary <- function() {
  ii <- pammtools::int_info(CUT)
  bind_rows(
    mutate(ii, x_num = 0, x = factor(0, levels = c(0, 1))),
    mutate(ii, x_num = 1, x = factor(1, levels = c(0, 1)))
  ) |>
    group_by(x) |>
    arrange(tend, .by_group = TRUE)
}

grid_x0 <- function() {
  pammtools::int_info(CUT) |>
    mutate(x = 0) |>
    group_by(x) |>
    arrange(tend, .by_group = TRUE)
}

grid_cr <- function() {
  ii <- pammtools::int_info(CUT)
  bind_rows(
    mutate(ii, x = 0, cause = factor(1, levels = c(1, 2))),
    mutate(ii, x = 0, cause = factor(2, levels = c(1, 2)))
  ) |>
    group_by(cause) |>
    arrange(tend, .by_group = TRUE)
}

make_grid <- function(cell) {
  switch(
    cell$model_type,
    ph = ,
    tv = grid_binary(),
    smooth = grid_x0(),
    cr = grid_cr()
  )
}

# ---- shared extraction helpers -------------------------------------------------
surv_rows_from <- function(sv, binary) {
  sv <- ungroup(sv) |> filter(tend %in% TEVAL)
  qoi_rows(
    "surv",
    sv$tend,
    if (binary) sv$x_num else sv$x,
    sv$surv_prob,
    (sv$surv_upper - sv$surv_lower) / (2 * Z975),
    sv$surv_lower,
    sv$surv_upper
  )
}

hazard_rows_from <- function(hz, binary) {
  hz <- ungroup(hz) |> filter(tend %in% TEVAL)
  qoi_rows(
    "hazard",
    hz$tend,
    if (binary) hz$x_num else hz$x,
    hz$hazard,
    (hz$ci_upper - hz$ci_lower) / (2 * Z975),
    hz$ci_lower,
    hz$ci_upper
  )
}

# lpmatrix difference (x=1 minus x=0) on the TEVAL intervals: logHR(t) = D beta
loghr_design <- function(fit) {
  ii <- pammtools::int_info(CUT) |> filter(tend %in% TEVAL)
  nd1 <- mutate(ii, x_num = 1, x = factor(1, levels = c(0, 1)))
  nd0 <- mutate(ii, x_num = 0, x = factor(0, levels = c(0, 1)))
  predict(fit, newdata = nd1, type = "lpmatrix") -
    predict(fit, newdata = nd0, type = "lpmatrix")
}

# logHR draws for a single gam fit
loghr_rows_gam <- function(fit, nsim) {
  D <- loghr_design(fit)
  B <- mvtnorm::rmvnorm(nsim, mean = coef(fit), sigma = fit[["Vp"]])
  M <- D %*% t(B)
  est <- drop(D %*% coef(fit))
  qoi_rows(
    "logHR",
    TEVAL,
    NA_real_,
    est,
    apply(M, 1, stats::sd),
    apply(M, 1, stats::quantile, ALPHA / 2),
    apply(M, 1, stats::quantile, 1 - ALPHA / 2)
  )
}

# pooled logHR draws for a pamm_ic fit: paired draws (X1 - X0) %*% beta_draw
# from the SAME coefficient draw, per-fit design matrices, Rubin-inflated on
# the logHR scale (mirrors pammtools' ic_ci_draws + rubin_inflate machinery)
loghr_rows_mi <- function(object, nsim) {
  fits <- object[["fits"]]
  m <- length(fits)
  per <- ceiling(nsim / m)
  pieces <- lapply(fits, function(f) {
    D <- loghr_design(f)
    B <- mvtnorm::rmvnorm(per, mean = coef(f), sigma = f[["Vp"]])
    list(draws = D %*% t(B), estimate = drop(D %*% coef(f)))
  })
  M <- do.call(
    cbind,
    rubin_inflate_local(
      lapply(pieces, `[[`, "draws"),
      lapply(pieces, `[[`, "estimate")
    )
  )
  est <- rowMeans(do.call(cbind, lapply(pieces, `[[`, "estimate")))
  qoi_rows(
    "logHR",
    TEVAL,
    NA_real_,
    est,
    apply(M, 1, stats::sd),
    apply(M, 1, stats::quantile, ALPHA / 2),
    apply(M, 1, stats::quantile, 1 - ALPHA / 2)
  )
}

# centered smooth f(x) on X_GRID for a single gam fit (or one imputation fit)
fx_design <- function(fit) {
  nd <- tibble(tend = CUT[2], x = X_GRID)
  X <- predict(fit, newdata = nd, type = "lpmatrix")
  keep <- grepl("^s\\(x\\)", colnames(X))
  X[, !keep] <- 0
  # center estimate over the evaluation grid (truth is centered the same way)
  sweep(X, 2, colMeans(X))
}

fx_rows_gam <- function(fit, nsim) {
  Xc <- fx_design(fit)
  B <- mvtnorm::rmvnorm(nsim, mean = coef(fit), sigma = fit[["Vp"]])
  M <- Xc %*% t(B)
  est <- drop(Xc %*% coef(fit))
  qoi_rows(
    "f_x",
    NA_real_,
    X_GRID,
    est,
    apply(M, 1, stats::sd),
    apply(M, 1, stats::quantile, ALPHA / 2),
    apply(M, 1, stats::quantile, 1 - ALPHA / 2)
  )
}

fx_rows_mi <- function(object, nsim) {
  fits <- object[["fits"]]
  per <- ceiling(nsim / length(fits))
  pieces <- lapply(fits, function(f) {
    Xc <- fx_design(f)
    B <- mvtnorm::rmvnorm(per, mean = coef(f), sigma = f[["Vp"]])
    list(draws = Xc %*% t(B), estimate = drop(Xc %*% coef(f)))
  })
  M <- do.call(
    cbind,
    rubin_inflate_local(
      lapply(pieces, `[[`, "draws"),
      lapply(pieces, `[[`, "estimate")
    )
  )
  est <- rowMeans(do.call(cbind, lapply(pieces, `[[`, "estimate")))
  qoi_rows(
    "f_x",
    NA_real_,
    X_GRID,
    est,
    apply(M, 1, stats::sd),
    apply(M, 1, stats::quantile, ALPHA / 2),
    apply(M, 1, stats::quantile, 1 - ALPHA / 2)
  )
}

beta_rows_gam <- function(fit, coef_name, estimand = "beta") {
  est <- unname(coef(fit)[coef_name])
  se <- sqrt(stats::vcov(fit)[coef_name, coef_name])
  qoi_rows(
    estimand,
    NA_real_,
    NA_real_,
    est,
    se,
    est - Z975 * se,
    est + Z975 * se
  )
}

beta_rows_mi <- function(object, coef_name, estimand = "beta") {
  est <- vapply(object$fits, function(f) unname(coef(f)[coef_name]), 0)
  vv <- vapply(object$fits, function(f) stats::vcov(f)[coef_name, coef_name], 0)
  nu_com <- max(object$n_obs - sum(object$fits[[1]]$edf), 1)
  rb <- rubin_scalar_br(est, vv, nu_com)
  qoi_rows(estimand, NA_real_, NA_real_, rb$est, rb$se, rb$lower, rb$upper)
}

cif_rows_from <- function(cf) {
  cf <- ungroup(cf) |> filter(tend %in% TEVAL)
  qoi_rows(
    paste0("cif", cf$cause),
    cf$tend,
    0,
    cf$cif,
    (cf$cif_upper - cf$cif_lower) / (2 * Z975),
    cf$cif_lower,
    cf$cif_upper
  )
}

# midpoint event time / status from (L, R] bounds (as in the prototype)
mid_time <- function(L, R) {
  ifelse(
    is.infinite(R),
    L,
    ifelse(L == R, L, pmin((L + pmin(R, HORIZON)) / 2, HORIZON))
  )
}
ev_from_R <- function(R) ifelse(is.infinite(R), 0L, 1L)

# ---- internal methods: MI / midpoint / oracle ----------------------------------
ic_formula_for <- function(cell) {
  switch(
    cell$model_type,
    ph = ,
    tv = Surv(L, R, type = "interval2") ~ x + x_num,
    smooth = Surv(L, R, type = "interval2") ~ x,
    cr = Surv(L, R, type = "interval2") ~ x
  )
}

extract_gam_single <- function(fit, cell) {
  binary <- cell$model_type %in% c("ph", "tv")
  grid <- make_grid(cell)
  sv <- pammtools::add_surv_prob(grid, fit, ci_type = "sim", nsim = NSIM_RUN)
  hz <- pammtools::add_hazard(grid, fit, ci_type = "sim", nsim = NSIM_RUN)
  out <- bind_rows(surv_rows_from(sv, binary), hazard_rows_from(hz, binary))
  if (binary) {
    out <- bind_rows(out, loghr_rows_gam(fit, NSIM_RUN))
  }
  if (cell$model_type == "ph") {
    out <- bind_rows(out, beta_rows_gam(fit, "x1"))
  }
  if (cell$model_type == "smooth") {
    out <- bind_rows(out, fx_rows_gam(fit, NSIM_RUN))
  }
  out
}

fit_mi <- function(dat, cell) {
  m <- if (!is.na(M_RUN)) M_RUN else cell$m
  if (cell$model_type == "cr") {
    fit <- pammtools::pamm_ic_cr(
      ic_formula_for(cell),
      dat$icd,
      cause = "cause",
      model_formula = model_formula_for("cr"),
      cut = CUT,
      m = m
    )
    return(bind_rows(
      cif_rows_from(pammtools::add_cif(grid_cr(), fit, nsim = NSIM_RUN)),
      beta_rows_mi(fit, "cause1:x", "beta1"),
      beta_rows_mi(fit, "cause2:x", "beta2")
    ))
  }
  fit <- pammtools::pamm_ic(
    ic_formula_for(cell),
    dat$icd,
    model_formula = model_formula_for(cell$model_type),
    cut = CUT,
    m = m
  )
  binary <- cell$model_type %in% c("ph", "tv")
  grid <- make_grid(cell)
  sv <- pammtools::add_surv_prob(grid, fit, nsim = NSIM_RUN)
  hz <- pammtools::add_hazard(grid, fit, nsim = NSIM_RUN)
  out <- bind_rows(surv_rows_from(sv, binary), hazard_rows_from(hz, binary))
  if (binary) {
    out <- bind_rows(out, loghr_rows_mi(fit, NSIM_RUN))
  }
  if (cell$model_type == "ph") {
    out <- bind_rows(out, beta_rows_mi(fit, "x1"))
  }
  if (cell$model_type == "smooth") {
    out <- bind_rows(out, fx_rows_mi(fit, NSIM_RUN))
  }
  out
}

ped_formula_for <- function(cell, time_var, status_var) {
  rhs <- switch(
    cell$model_type,
    ph = ,
    tv = "x + x_num",
    smooth = "x",
    cr = "x"
  )
  stats::as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, status_var, rhs))
}

fit_midpoint <- function(dat, cell) {
  icd <- dat$icd
  icd$tm <- mid_time(icd$L, icd$R)
  icd$ev <- ev_from_R(icd$R)
  if (cell$model_type == "cr") {
    # naive midpoint with single-imputed causes: unknown -> most frequent cause
    cs <- icd$cause
    if (anyNA(cs)) {
      obs <- cs[!is.na(cs) & cs > 0]
      cs[is.na(cs)] <- as.integer(names(which.max(table(obs))))
    }
    icd$.s <- ifelse(icd$ev == 0, 0L, cs)
    ped <- pammtools::as_ped(
      icd[icd$tm > 0, ],
      ped_formula_for(cell, "tm", ".s"),
      cut = CUT
    )
    fit <- pammtools::pamm(model_formula_for("cr"), data = ped)
    return(bind_rows(
      cif_rows_from(pammtools::add_cif(grid_cr(), fit, nsim = NSIM_RUN)),
      beta_rows_gam(fit, "cause1:x", "beta1"),
      beta_rows_gam(fit, "cause2:x", "beta2")
    ))
  }
  ped <- pammtools::as_ped(
    icd[icd$tm > 0, ],
    ped_formula_for(cell, "tm", "ev"),
    cut = CUT
  )
  fit <- pammtools::pamm(model_formula_for(cell$model_type), data = ped)
  extract_gam_single(fit, cell)
}

fit_oracle <- function(dat, cell) {
  sdf <- dat$exact
  if (cell$model_type == "cr") {
    sdf$.s <- sdf$cause # true causes (pre-masking)
    ped <- pammtools::as_ped(
      sdf,
      ped_formula_for(cell, "time", ".s"),
      cut = CUT
    )
    fit <- pammtools::pamm(model_formula_for("cr"), data = ped)
    return(bind_rows(
      cif_rows_from(pammtools::add_cif(grid_cr(), fit, nsim = NSIM_RUN)),
      beta_rows_gam(fit, "cause1:x", "beta1"),
      beta_rows_gam(fit, "cause2:x", "beta2")
    ))
  }
  ped <- pammtools::as_ped(
    sdf,
    ped_formula_for(cell, "time", "status"),
    cut = CUT
  )
  fit <- pammtools::pamm(model_formula_for(cell$model_type), data = ped)
  extract_gam_single(fit, cell)
}

# ---- external methods -----------------------------------------------------------
# icenReg Weibull PH: H(t|x) = (t/scale)^shape * exp(beta x)
fit_ic_par <- function(dat, cell) {
  fit <- icenReg::ic_par(
    Surv(L, R, type = "interval2") ~ x_num,
    data = dat$icd,
    model = "ph",
    dist = "weibull"
  )
  cf <- fit$coefficients
  # qoi_par() indexes positionally: assert ORDER, not just membership
  stopifnot(identical(names(cf)[1:3], c("log_shape", "log_scale", "x_num")))
  V <- fit$var
  ti <- teval_intervals()
  # icenReg centers covariates internally: baseline = hazard at the covariate
  # MEAN; un-center via covarOffset (verified against getFitEsts in P1)
  off <- as.numeric(fit$covarOffset)

  qoi_par <- function(b) {
    # b = c(log_shape, log_scale, beta); H0(t) = (t/scale)^shape
    H <- function(t, xv) {
      ifelse(t <= 0, 0, exp(exp(b[1]) * (log(pmax(t, 1e-12)) - b[2]))) *
        exp(b[3] * (xv - off))
    }
    c(
      surv = c(exp(-H(TEVAL, 0)), exp(-H(TEVAL, 1))),
      hazard = c(
        (H(ti$tend, 0) - H(ti$tstart, 0)) / ti$intlen,
        (H(ti$tend, 1) - H(ti$tstart, 1)) / ti$intlen
      ),
      loghr = rep(b[3], length(TEVAL))
    )
  }
  est <- qoi_par(cf)
  B <- mvtnorm::rmvnorm(NSIM_RUN, mean = cf, sigma = V)
  M <- apply(B, 1, qoi_par) # n_qoi x NSIM
  lo <- apply(M, 1, stats::quantile, ALPHA / 2)
  hi <- apply(M, 1, stats::quantile, 1 - ALPHA / 2)
  sd_ <- apply(M, 1, stats::sd)

  k <- length(TEVAL)
  i_sv <- seq_len(2 * k)
  i_hz <- 2 * k + seq_len(2 * k)
  i_lh <- 4 * k + seq_len(k)
  se_b <- sqrt(V["x_num", "x_num"])
  b <- unname(cf["x_num"])
  out <- bind_rows(
    qoi_rows(
      "surv",
      rep(TEVAL, 2),
      rep(c(0, 1), each = k),
      est[i_sv],
      sd_[i_sv],
      lo[i_sv],
      hi[i_sv]
    ),
    qoi_rows(
      "hazard",
      rep(TEVAL, 2),
      rep(c(0, 1), each = k),
      est[i_hz],
      sd_[i_hz],
      lo[i_hz],
      hi[i_hz]
    ),
    qoi_rows("logHR", TEVAL, NA_real_, est[i_lh], sd_[i_lh], lo[i_lh], hi[i_lh])
  )
  if (cell$model_type == "ph") {
    out <- bind_rows(
      out,
      qoi_rows(
        "beta",
        NA_real_,
        NA_real_,
        b,
        se_b,
        b - Z975 * se_b,
        b + Z975 * se_b
      )
    )
  }
  out
}

# icenReg semiparametric Cox PH: beta everywhere; bootstrap CI in subset cells;
# S(t|x) point estimates only (NPMLE baseline has no plug-in CI)
fit_ic_sp <- function(dat, cell) {
  bs <- if (isTRUE(cell$icsp_boot)) BS_SAMPLES_ICSP else 0L
  fit <- icenReg::ic_sp(
    Surv(L, R, type = "interval2") ~ x_num,
    data = dat$icd,
    bs_samples = bs
  )
  b <- unname(fit$coefficients["x_num"])
  if (bs > 0) {
    bsamp <- fit$bsMat[, "x_num"]
    bsamp <- bsamp[is.finite(bsamp)] # failed bootstrap draws may be NA
    blo <- unname(stats::quantile(bsamp, ALPHA / 2))
    bhi <- unname(stats::quantile(bsamp, 1 - ALPHA / 2))
    bse <- stats::sd(bsamp)
  } else {
    blo <- bhi <- bse <- NA_real_
  }
  s0 <- 1 - icenReg::getFitEsts(fit, newdata = data.frame(x_num = 0), q = TEVAL)
  s1 <- 1 - icenReg::getFitEsts(fit, newdata = data.frame(x_num = 1), q = TEVAL)
  k <- length(TEVAL)
  out <- bind_rows(
    qoi_rows(
      "surv",
      rep(TEVAL, 2),
      rep(c(0, 1), each = k),
      c(s0, s1),
      NA_real_,
      NA_real_,
      NA_real_
    ),
    qoi_rows("logHR", TEVAL, NA_real_, rep(b, k), bse, blo, bhi)
  )
  if (cell$model_type == "ph") {
    out <- bind_rows(
      out,
      qoi_rows("beta", NA_real_, NA_real_, b, bse, blo, bhi)
    )
  }
  out
}

# Turnbull NPMLE per stratum; point accuracy only. The NPMLE is defined only up
# to equivalence inside Turnbull intervals: score the midpoint of the upper and
# lower step curves, record the half-gap as `ambiguity`.
turnbull_S_at <- function(fit_np, times) {
  sc <- icenReg::getSCurves(fit_np)
  tb <- sc$Tbull_ints # matrix: cols = lower/upper interval endpoints
  s <- sc$S_curves[[1]]
  # upper envelope: S only drops at the RIGHT endpoint of each Turnbull
  # interval (holds the pre-drop value through it); lower envelope drops
  # already AT the left endpoint
  up_fun <- stats::stepfun(tb[, 2], c(1, s), right = TRUE)
  lo_fun <- stats::stepfun(tb[, 1], c(1, s), right = FALSE)
  up <- up_fun(times)
  loval <- lo_fun(times)
  stopifnot(all(up - loval >= -1e-12))
  tibble(est = (up + loval) / 2, ambiguity = (up - loval) / 2)
}

fit_turnbull <- function(dat, cell) {
  icd <- dat$icd
  res <- lapply(c(0, 1), function(xv) {
    sub <- icd[icd$x_num == xv, c("L", "R")]
    fit <- icenReg::ic_np(cbind(L, R) ~ 0, data = sub)
    s <- turnbull_S_at(fit, TEVAL)
    qoi_rows(
      "surv",
      TEVAL,
      xv,
      s$est,
      NA_real_,
      NA_real_,
      NA_real_,
      ambiguity = s$ambiguity
    )
  })
  bind_rows(res)
}

# ---- dispatcher --------------------------------------------------------------------
FIT_FUNS <- list(
  mi = fit_mi,
  midpoint = fit_midpoint,
  oracle = fit_oracle,
  ic_par = fit_ic_par,
  ic_sp = fit_ic_sp,
  turnbull = fit_turnbull
)

# misspecified-by-design flags (per method x cell x estimand)
flag_misspecified <- function(method, cell, estimand) {
  if (method == "ic_par") {
    bad_baseline <- cell$baseline != "const"
    tv <- cell$effect == "tv"
    return(
      (estimand %in% c("surv", "hazard") & (bad_baseline | tv)) |
        (estimand %in% c("logHR", "beta") & (tv | bad_baseline))
    )
  }
  if (method == "ic_sp") {
    tv <- cell$effect == "tv"
    return(estimand %in% c("logHR", "beta", "surv") & tv)
  }
  rep(FALSE, length(estimand))
}

# Fit one method with per-method sub-seed, warning capture and timing.
# Never errors: failures yield a single row with error_msg set.
run_method <- function(method, dat, cell, subseed) {
  set.seed(subseed)
  n_warn <- 0L
  first_warn <- NA_character_
  t0 <- proc.time()[3]
  est <- tryCatch(
    withCallingHandlers(
      FIT_FUNS[[method]](dat, cell),
      warning = function(w) {
        n_warn <<- n_warn + 1L
        if (is.na(first_warn)) first_warn <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      tibble(
        estimand = NA_character_,
        t = NA_real_,
        x = NA_real_,
        est = NA_real_,
        se = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        ambiguity = NA_real_,
        error_msg = conditionMessage(e)
      )
    }
  )
  if (!"error_msg" %in% names(est)) {
    est$error_msg <- NA_character_
  }
  est$method <- method
  est$fit_time <- proc.time()[3] - t0
  est$n_warnings <- n_warn
  est$first_warning <- first_warn
  est$misspecified <- flag_misspecified(method, cell, est$estimand)
  est
}
