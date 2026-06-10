# ===========================================================================
# IC benchmark: data-generating processes and exact truths
# ---------------------------------------------------------------------------
# sim_pexp() draws from a LEFT-ENDPOINT piecewise-constant hazard on FINE.
# All truths here are computed exactly from those same rate vectors
# (cumulative sums; H(t) is piecewise linear), so truth and simulator are
# derived from one definition and cannot diverge.
# Requires config.R to be sourced first.
# ===========================================================================

# ---- log-hazard components --------------------------------------------------
lb_const <- function(t) rep(-1.5, length(t))
lb_peaked <- function(t) -3.2 + 1.6 * exp(-0.5 * ((t - 3) / 1)^2)
lb_rising <- function(t) -3.7 + 0.12 * t
beta_tv <- function(t) 1.2 * exp(-0.3 * t)
f_smooth <- function(x) 0.8 * sin(pi * x)
# competing risks (prototype truths)
lb_cr1 <- function(t) -3.3 + 1.2 * exp(-0.5 * ((t - 2.5) / 1)^2)
lb_cr2 <- function(t) -3.7 + 0.12 * t

baseline_fun <- function(baseline) {
  switch(
    baseline,
    const = lb_const,
    peaked = lb_peaked,
    rising = lb_rising,
    stop("unknown baseline: ", baseline)
  )
}

# log-hazard as a function of (t, x) for single-event cells
cell_loghaz <- function(cell) {
  lb <- baseline_fun(cell$baseline)
  switch(
    cell$effect,
    ph = function(t, x) lb(t) + BETA_PH * x,
    tv = function(t, x) lb(t) + beta_tv(t) * x,
    smooth = function(t, x) lb(t) + f_smooth(x),
    stop("unknown effect: ", cell$effect)
  )
}

# sim_pexp formulas built from the SAME functional forms (sim_pexp evaluates
# these at the left endpoints of FINE, i.e. exactly the rate vectors below)
cell_sim_formula <- function(cell) {
  base_str <- switch(
    cell$baseline,
    const = "-1.5",
    peaked = "-3.2 + 1.6 * exp(-0.5 * ((t - 3) / 1)^2)",
    rising = "-3.7 + 0.12 * t"
  )
  eff_str <- switch(
    cell$effect,
    ph = sprintf("%s * x_num", format(BETA_PH)),
    tv = "(1.2 * exp(-0.3 * t)) * x_num",
    smooth = "0.8 * sin(pi * x)"
  )
  stats::as.formula(paste("~", base_str, "+", eff_str))
}

CR_SIM_F1 <- ~ -3.3 + 1.2 * exp(-0.5 * ((t - 2.5) / 1)^2) + 0.7 * x
CR_SIM_F2 <- ~ -3.7 + 0.12 * t - 0.5 * x

# ---- exact cumulative hazard machinery ---------------------------------------
# rates on FINE intervals (left-endpoint evaluation, as in sim_pexp)
fine_rates <- function(loghaz, x) {
  left <- FINE[-length(FINE)]
  exp(loghaz(left, x))
}

# exact H(t) for a piecewise-constant rate vector: piecewise-linear interpolation
cumhaz_fun <- function(rates) {
  H <- c(0, cumsum(rates * FINE_STEP))
  stats::approxfun(FINE, H, rule = 2)
}

# ---- truth tables -------------------------------------------------------------
# CUT intervals whose tend is in TEVAL (hazard/logHR are interval-average there)
teval_intervals <- function() {
  j <- match(TEVAL, CUT)
  tibble(tstart = CUT[j - 1], tend = CUT[j], intlen = CUT[j] - CUT[j - 1])
}

# Truth rows for a single-event binary-x cell:
#   surv:   S(t | x) at TEVAL x {0,1}
#   hazard: interval-average  hbar_j(x) = (H(tend) - H(tstart)) / intlen
#   logHR:  log hbar_j(1) - log hbar_j(0)   (constant = BETA_PH under PH)
#   beta:   scalar (PH cells only)
truth_single <- function(cell) {
  lh <- cell_loghaz(cell)
  H0 <- cumhaz_fun(fine_rates(lh, 0))
  H1 <- cumhaz_fun(fine_rates(lh, 1))
  ti <- teval_intervals()
  hbar0 <- (H0(ti$tend) - H0(ti$tstart)) / ti$intlen
  hbar1 <- (H1(ti$tend) - H1(ti$tstart)) / ti$intlen

  out <- bind_rows(
    tibble(estimand = "surv", t = TEVAL, x = 0, truth = exp(-H0(TEVAL))),
    tibble(estimand = "surv", t = TEVAL, x = 1, truth = exp(-H1(TEVAL))),
    tibble(estimand = "hazard", t = TEVAL, x = 0, truth = hbar0),
    tibble(estimand = "hazard", t = TEVAL, x = 1, truth = hbar1),
    tibble(
      estimand = "logHR",
      t = TEVAL,
      x = NA_real_,
      truth = log(hbar1) - log(hbar0)
    )
  )
  if (cell$effect == "ph") {
    out <- bind_rows(
      out,
      tibble(estimand = "beta", t = NA_real_, x = NA_real_, truth = BETA_PH)
    )
  }
  out
}

# Spoke S: surv/hazard at x = 0 plus centered f(x) on X_GRID
truth_smooth <- function(cell) {
  lh <- cell_loghaz(cell)
  H0 <- cumhaz_fun(fine_rates(lh, 0))
  ti <- teval_intervals()
  fx <- f_smooth(X_GRID)
  bind_rows(
    tibble(estimand = "surv", t = TEVAL, x = 0, truth = exp(-H0(TEVAL))),
    tibble(
      estimand = "hazard",
      t = TEVAL,
      x = 0,
      truth = (H0(ti$tend) - H0(ti$tstart)) / ti$intlen
    ),
    tibble(estimand = "f_x", t = NA_real_, x = X_GRID, truth = fx - mean(fx)) # estimate is centered the same way
  )
}

# CR arm: cause-specific betas + CIF_k(t) at x = 0, computed exactly from the
# piecewise-constant cause-specific rates:
#   CIF_k(t) = sum_j  r_kj / r_.j * S_all(t_{j-1}) * (1 - exp(-r_.j * dt_j))
truth_cr <- function(cell) {
  r1 <- fine_rates(function(t, x) lb_cr1(t) + B1_CR * x, 0)
  r2 <- fine_rates(function(t, x) lb_cr2(t) + B2_CR * x, 0)
  rtot <- r1 + r2
  Sall_left <- exp(-c(0, cumsum(rtot * FINE_STEP)))[-(length(FINE))]
  inc <- Sall_left * (1 - exp(-rtot * FINE_STEP))
  cif1 <- c(0, cumsum((r1 / rtot) * inc))
  cif2 <- c(0, cumsum((r2 / rtot) * inc))
  cif1_at <- stats::approx(FINE, cif1, TEVAL)$y # exact at TEVAL (in FINE)
  cif2_at <- stats::approx(FINE, cif2, TEVAL)$y
  bind_rows(
    tibble(estimand = "beta1", t = NA_real_, x = NA_real_, truth = B1_CR),
    tibble(estimand = "beta2", t = NA_real_, x = NA_real_, truth = B2_CR),
    tibble(estimand = "cif1", t = TEVAL, x = 0, truth = cif1_at),
    tibble(estimand = "cif2", t = TEVAL, x = 0, truth = cif2_at)
  )
}

make_truth_table <- function(cell) {
  switch(
    cell$model_type,
    ph = ,
    tv = truth_single(cell),
    smooth = truth_smooth(cell),
    cr = truth_cr(cell)
  )
}

# beta(t) reference curve for TV cells (reporting / projection error)
true_beta_t <- beta_tv

# ---- data generation -----------------------------------------------------------
# Returns list(icd = interval-censored data, exact = exact-time data, meta = tibble)
generate_data <- function(cell, seed_data) {
  set.seed(seed_data)
  # reserve the first draws for method sub-seeds (must match method_subseeds())
  invisible(sample.int(.Machine$integer.max, 32L))

  n <- cell$n
  if (cell$model_type %in% c("ph", "tv")) {
    d <- tibble(
      id = seq_len(n),
      x_num = stats::rbinom(n, 1, 0.5)
    ) |>
      mutate(x = factor(x_num, levels = c(0, 1)))
    sdf <- pammtools::sim_pexp(cell_sim_formula(cell), d, cut = FINE)
  } else if (cell$model_type == "smooth") {
    d <- tibble(id = seq_len(n), x = stats::runif(n, -1, 1))
    sdf <- pammtools::sim_pexp(cell_sim_formula(cell), d, cut = FINE)
  } else {
    # cr: two independent latent cause-specific times
    d <- tibble(id = seq_len(n), x = stats::runif(n, -1, 1))
    T1 <- pammtools::sim_pexp(CR_SIM_F1, d, cut = FINE)$time
    T2 <- pammtools::sim_pexp(CR_SIM_F2, d, cut = FINE)$time
    obs <- pmin(T1, T2, HORIZON)
    cause <- ifelse(obs >= HORIZON, 0L, ifelse(T1 < T2, 1L, 2L))
    sdf <- d |>
      mutate(time = obs, status = as.integer(cause > 0), cause = cause)
  }

  icd <- pammtools::add_inspections(
    sdf,
    mechanism = cell$mechanism,
    rate = cell$rate,
    schedule = cell$schedule[[1]],
    max_time = HORIZON
  )
  if (cell$model_type == "cr" && cell$cause_mask > 0) {
    ev <- which(icd$cause > 0)
    mask <- ev[sample.int(length(ev), round(cell$cause_mask * length(ev)))]
    icd$cause[mask] <- NA_integer_
  }

  fin <- is.finite(icd$R) & icd$R > icd$L
  meta <- tibble(
    mean_int = mean(pmin(icd$R, HORIZON)[fin] - icd$L[fin]),
    prop_right_cens = mean(is.infinite(icd$R)),
    prop_left_cens = mean(icd$L == 0 & is.finite(icd$R))
  )
  list(icd = icd, exact = sdf, meta = meta)
}

# ---- least-false constant logHR for TV cells (reporting overlay) ----------------
# Projection of the TV truth onto a constant-beta Cox / Weibull PH model, via one
# large-n fit on exact times admin-censored at HORIZON. Convention: exact-time
# projection (cheap, censoring-scheme-free); documented in README.
compute_least_false <- function(cell, n = 2e5, seed = 1L) {
  stopifnot(cell$effect == "tv")
  set.seed(seed)
  d <- tibble(id = seq_len(n), x_num = stats::rbinom(n, 1, 0.5)) |>
    mutate(x = factor(x_num, levels = c(0, 1)))
  sdf <- pammtools::sim_pexp(cell_sim_formula(cell), d, cut = FINE)
  cox <- survival::coxph(survival::Surv(time, status) ~ x_num, data = sdf)
  wei <- survival::survreg(
    survival::Surv(time, status) ~ x_num,
    data = sdf,
    dist = "weibull"
  )
  # survreg AFT -> PH: beta_PH = -beta_AFT / scale (valid for Weibull)
  tibble(
    cell_id = cell$cell_id,
    lf_cox = unname(stats::coef(cox)["x_num"]),
    lf_weibull = unname(-stats::coef(wei)["x_num"] / wei$scale)
  )
}
