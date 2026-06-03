#!/usr/bin/env Rscript
# ===========================================================================
# Interval-censored PAMMs: simulation study
# ---------------------------------------------------------------------------
# Demonstrates, on clean synthetic data, for BOTH single-event and competing
# risks, that the multiple-imputation (MI) workflow in pammtools
# (pamm_ic / pamm_ic_cr) improves on naive single midpoint imputation:
#
#   (A) CONFIDENCE-INTERVAL COVERAGE of the time-dependent functionals that
#       practitioners report -- survival S(t) and cumulative incidence CIF_k(t)
#       -- is restored to ~nominal by MI, whereas naive midpoint imputation
#       badly under-covers (it ignores the uncertainty about the unobserved
#       exact event times), increasingly so as inspection intervals lengthen.
#       The covariate-effect (slope) CIs, by contrast, are ~nominal under BOTH
#       methods: the regression slope is robust to interval censoring; the
#       information loss falls on the time dimension (baseline / S(t) / CIF).
#   (B) POINT ESTIMATES are no worse than midpoint and improve under heavy
#       interval censoring (long inspection intervals), approaching the
#       (unattainable) "oracle" fit on the exact event times.
#
# An oracle fit on the exact times is included throughout as the gold-standard
# lower bound. Results are averaged over many replications across a gradient of
# inspection densities (mean inspection-interval length).
#
# Run (from the package root):
#   IC_NREP=200 Rscript attic/ic-simulation-study.R
# A smaller/faster run:  IC_NREP=40 Rscript attic/ic-simulation-study.R
# Figures are written to attic/figures/.
#
# NOTE: this is a stand-alone study script kept in attic/, not part of the
# package or its test suite. Requires the package to be installed/loadable plus
# ggplot2, dplyr, tidyr, purrr (and parallel, used if >1 core is available).
# ===========================================================================

suppressMessages({
  library(pammtools)
  library(survival)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

# ---- configuration --------------------------------------------------------
NREP    <- as.integer(Sys.getenv("IC_NREP", "200"))
NCORES  <- as.integer(Sys.getenv("IC_NCORES",
  as.character(max(1L, parallel::detectCores() - 1L))))
M       <- 10L                       # imputations per MI fit
NSIM    <- 200L                      # posterior draws for functional CIs
HORIZON <- 10
FINE    <- seq(0, HORIZON, by = 0.05)
CUT     <- seq(0, HORIZON, by = 0.5) # analysis interval grid (fixed, shared)
TEVAL   <- c(1, 2, 3, 5, 7)          # times at which S(t)/CIF are evaluated
RATES   <- c(1.5, 0.6, 0.3)          # inspection rates: dense -> sparse
FIGDIR  <- file.path("attic", "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(20240603)
theme_set(theme_bw(base_size = 12))
cat(sprintf("Running IC simulation study: NREP=%d, cores=%d\n", NREP, NCORES))

# ---- truth: single event --------------------------------------------------
# log-hazard = peaked baseline + linear covariate effect (truth beta_x)
BETA_X    <- 0.6
lbase_se  <- function(t) -3.2 + 1.6 * exp(-0.5 * ((t - 3) / 1.0)^2)
SIM_SE    <- ~ -3.2 + 1.6 * exp(-0.5 * ((t - 3) / 1.0)^2) + 0.6 * x
haz_se    <- function(t, x = 0) exp(lbase_se(t) + BETA_X * x)
S_se      <- function(t, x = 0)
  exp(-vapply(t, function(tt) integrate(haz_se, 0, tt, x = x)$value, numeric(1)))
TRUE_S    <- S_se(TEVAL, 0)

# ---- truth: competing risks (two causes, independent latent times) --------
B1 <- 0.7; B2 <- -0.5
lb1 <- function(t) -3.3 + 1.2 * exp(-0.5 * ((t - 2.5) / 1.0)^2)  # cause 1: peaked
lb2 <- function(t) -3.7 + 0.12 * t                               # cause 2: rising
SIM_C1 <- ~ -3.3 + 1.2 * exp(-0.5 * ((t - 2.5) / 1.0)^2) + 0.7 * x
SIM_C2 <- ~ -3.7 + 0.12 * t - 0.5 * x
h1 <- function(t, x = 0) exp(lb1(t) + B1 * x)
h2 <- function(t, x = 0) exp(lb2(t) + B2 * x)
# CIF_k(t) = int_0^t h_k(u) S_all(u) du, S_all = exp(-int (h1+h2))
.fine_u  <- FINE
.Sall_u  <- exp(-cumsum(c(0, (h1(.fine_u, 0) + h2(.fine_u, 0))[-1] * diff(.fine_u))))
.cif_of  <- function(hf) cumsum(c(0, (hf(.fine_u, 0) * .Sall_u)[-1] * diff(.fine_u)))
.cif1_u  <- .cif_of(h1); .cif2_u <- .cif_of(h2)
TRUE_CIF1 <- approx(.fine_u, .cif1_u, TEVAL)$y
TRUE_CIF2 <- approx(.fine_u, .cif2_u, TEVAL)$y

# ---- helpers --------------------------------------------------------------
# Rubin's rules for a scalar quantity pooled over m imputations.
rubin_scalar <- function(est, var, level = 0.95) {
  m <- length(est); Qbar <- mean(est); W <- mean(var); B <- stats::var(est)
  Tot <- W + (1 + 1 / m) * B; se <- sqrt(Tot)
  riv <- (1 + 1 / m) * B / W
  df  <- if (riv <= 0) 1e6 else (m - 1) * (1 + 1 / riv)^2
  tc  <- stats::qt(1 - (1 - level) / 2, df)
  c(est = Qbar, lower = Qbar - tc * se, upper = Qbar + tc * se)
}
normal_ci <- function(est, se, level = 0.95) {
  est <- unname(est); se <- unname(se)
  z <- stats::qnorm(1 - (1 - level) / 2)
  c(est = est, lower = est - z * se, upper = est + z * se)
}
covered <- function(ci, truth) unname(truth >= ci["lower"] & truth <= ci["upper"])

# midpoint event time / status from (L, R] bounds (no package internals)
mid_time   <- function(L, R) ifelse(is.infinite(R), L,
  ifelse(L == R, L, pmin((L + pmin(R, HORIZON)) / 2, HORIZON)))
ev_from_R  <- function(R) ifelse(is.infinite(R), 0L, 1L)
mean_int   <- function(L, R) mean((pmin(R, HORIZON) - L)[is.finite(R) & R > L])

row_term <- function(task, rate, mint, method, quantity, truth, ci) {
  data.frame(task = task, rate = rate, mean_int = mint, method = method,
    quantity = quantity, time = NA_real_, truth = truth,
    est = unname(ci["est"]), covered = covered(ci, truth),
    width = unname(ci["upper"] - ci["lower"]),
    sqerr = (unname(ci["est"]) - truth)^2, stringsAsFactors = FALSE)
}
# prediction grids built from the fixed cut grid so they always span [0, HORIZON]
grid_single <- function() transform(int_info(CUT), x = 0)
grid_cr <- function() {
  ii <- int_info(CUT)
  g  <- do.call(rbind, lapply(c("1", "2"),
    function(cl) transform(ii, x = 0, cause = factor(cl, levels = c("1", "2")))))
  dplyr::group_by(g, cause)
}

row_fun <- function(task, rate, mint, method, quantity, df, truth, lo, hi, est) {
  data.frame(task = task, rate = rate, mean_int = mint, method = method,
    quantity = quantity, time = df, truth = truth,
    est = est, covered = truth >= lo & truth <= hi, width = hi - lo,
    sqerr = (est - truth)^2, stringsAsFactors = FALSE)
}

# =====================  single-event replicate  ============================
one_rep_single <- function(rate, seed, n = 300) {
  set.seed(seed)
  d   <- data.frame(id = seq_len(n), x = runif(n, -1, 1))
  sdf <- sim_pexp(SIM_SE, d, cut = FINE)
  icd <- add_inspections(sdf, rate = rate, max_time = HORIZON)
  mint <- mean_int(icd$L, icd$R)
  f <- Surv(L, R, type = "interval2") ~ x

  fit <- suppressWarnings(pamm_ic(f, icd, cut = CUT, m = M))
  icd$tm <- mid_time(icd$L, icd$R); icd$ev <- ev_from_R(icd$R)
  mid <- suppressWarnings(pamm(ped_status ~ s(tend) + x,
    data = as_ped(icd[icd$tm > 0, ], Surv(tm, ev) ~ x, cut = CUT)))
  ora <- suppressWarnings(pamm(ped_status ~ s(tend) + x,
    data = as_ped(sdf, Surv(time, status) ~ x, cut = CUT)))

  ## --- term: beta_x ---
  ci_mi  <- rubin_scalar(vapply(fit$fits, function(z) coef(z)["x"], 0),
    vapply(fit$fits, function(z) vcov(z)["x", "x"], 0))
  ci_md  <- normal_ci(coef(mid)["x"], sqrt(vcov(mid)["x", "x"]))
  ci_or  <- normal_ci(coef(ora)["x"], sqrt(vcov(ora)["x", "x"]))
  rows <- list(
    row_term("single", rate, mint, "MI",       "beta[x]", BETA_X, ci_mi),
    row_term("single", rate, mint, "midpoint", "beta[x]", BETA_X, ci_md),
    row_term("single", rate, mint, "oracle",   "beta[x]", BETA_X, ci_or))

  ## --- functional: S(t) ---
  nd <- grid_single()
  grab <- function(s, lo, hi) {
    i <- match(TEVAL, s$tend)
    list(est = s[[ "surv_prob" ]][i], lo = s[[lo]][i], hi = s[[hi]][i])
  }
  sMI <- add_surv_prob(nd, fit, nsim = NSIM)
  sMD <- add_surv_prob(nd, mid, ci_type = "sim", nsim = NSIM)
  sOR <- add_surv_prob(nd, ora, ci_type = "sim", nsim = NSIM)
  gMI <- grab(sMI, "surv_lower", "surv_upper")
  gMD <- grab(sMD, "surv_lower", "surv_upper")
  gOR <- grab(sOR, "surv_lower", "surv_upper")
  rows <- c(rows, list(
    row_fun("single", rate, mint, "MI",       "S(t)", TEVAL, TRUE_S, gMI$lo, gMI$hi, gMI$est),
    row_fun("single", rate, mint, "midpoint", "S(t)", TEVAL, TRUE_S, gMD$lo, gMD$hi, gMD$est),
    row_fun("single", rate, mint, "oracle",   "S(t)", TEVAL, TRUE_S, gOR$lo, gOR$hi, gOR$est)))
  do.call(rbind, rows)
}

# =====================  competing-risks replicate  =========================
CR_FORM <- ped_status ~ s(tend, by = cause) + cause + cause:x

cr_data <- function(seed, n) {
  set.seed(seed)
  d  <- data.frame(id = seq_len(n), x = runif(n, -1, 1))
  T1 <- sim_pexp(SIM_C1, d, cut = FINE)$time
  T2 <- sim_pexp(SIM_C2, d, cut = FINE)$time
  obs <- pmin(T1, T2, HORIZON)
  cause <- ifelse(obs >= HORIZON, 0L, ifelse(T1 < T2, 1L, 2L))
  data.frame(id = d$id, x = d$x, time = obs,
    status = as.integer(cause > 0), cause = cause)
}
one_rep_cr <- function(rate, seed, n = 450) {
  dat <- cr_data(seed, n)
  icd <- add_inspections(dat, rate = rate, max_time = HORIZON)  # uses time,status
  icd$cause <- dat$cause
  mint <- mean_int(icd$L, icd$R)
  f <- Surv(L, R, type = "interval2") ~ x

  fit <- suppressWarnings(pamm_ic_cr(f, icd, cause = "cause",
    model_formula = CR_FORM, cut = CUT, m = M))
  icd$tm <- mid_time(icd$L, icd$R); icd$ev <- ev_from_R(icd$R)
  md_dat <- icd[icd$tm > 0, ]
  md_status <- ifelse(md_dat$ev == 0, 0L, md_dat$cause)
  ped_md <- as_ped(transform(md_dat, .t = tm, .s = md_status),
    Surv(.t, .s) ~ x, cut = CUT)
  mid <- suppressWarnings(pamm(CR_FORM, data = ped_md))
  ped_or <- as_ped(transform(dat, .s = cause), Surv(time, .s) ~ x, cut = CUT)
  ora <- suppressWarnings(pamm(CR_FORM, data = ped_or))

  ## --- terms: cause-specific slopes ---
  cs <- function(fit_obj, nm) {
    e <- vapply(fit_obj$fits, function(z) coef(z)[nm], 0)
    v <- vapply(fit_obj$fits, function(z) vcov(z)[nm, nm], 0)
    rubin_scalar(e, v)
  }
  rows <- list(
    row_term("cr", rate, mint, "MI", "beta[1]", B1, cs(fit, "cause1:x")),
    row_term("cr", rate, mint, "MI", "beta[2]", B2, cs(fit, "cause2:x")),
    row_term("cr", rate, mint, "midpoint", "beta[1]", B1,
      normal_ci(coef(mid)["cause1:x"], sqrt(vcov(mid)["cause1:x", "cause1:x"]))),
    row_term("cr", rate, mint, "midpoint", "beta[2]", B2,
      normal_ci(coef(mid)["cause2:x"], sqrt(vcov(mid)["cause2:x", "cause2:x"]))),
    row_term("cr", rate, mint, "oracle", "beta[1]", B1,
      normal_ci(coef(ora)["cause1:x"], sqrt(vcov(ora)["cause1:x", "cause1:x"]))),
    row_term("cr", rate, mint, "oracle", "beta[2]", B2,
      normal_ci(coef(ora)["cause2:x"], sqrt(vcov(ora)["cause2:x", "cause2:x"]))))

  ## --- functional: CIF_k(t) ---
  nd <- grid_cr()
  cifrows <- function(cf, method) {
    out <- list()
    for (cl in c(1, 2)) {
      sub <- cf[as.character(cf$cause) == as.character(cl), ]
      i <- match(TEVAL, sub$tend)
      truth <- if (cl == 1) TRUE_CIF1 else TRUE_CIF2
      out[[length(out) + 1]] <- row_fun("cr", rate, mint, method,
        paste0("CIF[", cl, "]"), TEVAL, truth,
        sub$cif_lower[i], sub$cif_upper[i], sub$cif[i])
    }
    out
  }
  cfMI <- add_cif(nd, fit, nsim = NSIM)
  cfMD <- add_cif(nd, mid, nsim = NSIM)
  cfOR <- add_cif(nd, ora, nsim = NSIM)
  rows <- c(rows, cifrows(cfMI, "MI"), cifrows(cfMD, "midpoint"),
    cifrows(cfOR, "oracle"))
  do.call(rbind, rows)
}

# ---- run ------------------------------------------------------------------
run_grid <- function(fun, n, tag) {
  jobs <- expand.grid(rate = RATES, rep = seq_len(NREP))
  seeds <- as.integer(1e6 + seq_len(nrow(jobs)))
  do_one <- function(i) {
    out <- tryCatch(fun(jobs$rate[i], seeds[i], n = n), error = function(e) NULL)
    if (!is.null(out)) out
  }
  cat(sprintf("  [%s] %d jobs ...\n", tag, nrow(jobs)))
  res <- if (NCORES > 1L) {
    parallel::mclapply(seq_len(nrow(jobs)), do_one, mc.cores = NCORES,
      mc.preschedule = FALSE)
  } else lapply(seq_len(nrow(jobs)), do_one)
  bind_rows(res)
}

RAW <- file.path(FIGDIR, "ic-sim-raw.rds")
if (nzchar(Sys.getenv("IC_REUSE")) && file.exists(RAW)) {
  res <- readRDS(RAW)
  cat(sprintf("Reusing cached raw results (%d rows); regenerating figures only.\n",
    nrow(res)))
} else {
  t0 <- Sys.time()
  res_se <- run_grid(one_rep_single, 300, "single")
  res_cr <- run_grid(one_rep_cr,     450, "cr")
  res <- bind_rows(res_se, res_cr)
  saveRDS(res, RAW)
  cat(sprintf("Done in %.1f min; %d result rows.\n",
    as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(res)))
}

# ---- aggregate ------------------------------------------------------------
agg <- res %>%
  group_by(task, quantity, method, rate) %>%
  summarise(
    mean_int = mean(mean_int, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    rmse     = sqrt(mean(sqerr, na.rm = TRUE)),
    bias     = mean(est - truth, na.rm = TRUE),
    width    = mean(width, na.rm = TRUE),
    nrep     = dplyr::n(),
    .groups  = "drop") %>%
  mutate(method = factor(method, levels = c("oracle", "MI", "midpoint")))
write.csv(agg, file.path(FIGDIR, "ic-sim-summary.csv"), row.names = FALSE)
show_cov <- function(keep) print(as.data.frame(agg %>% filter(keep) %>%
  mutate(mean_int = round(mean_int, 2), coverage = round(coverage, 3)) %>%
  arrange(task, quantity, desc(mean_int), method) %>%
  select(task, quantity, mean_int, method, coverage)), row.names = FALSE)
cat("\n==== coverage of FUNCTIONALS S(t)/CIF (target 0.95) ====\n")
cat("    [headline: midpoint under-covers under heavy IC; MI tracks the oracle]\n")
show_cov(!grepl("beta", agg$quantity))
cat("\n==== coverage of covariate TERMS (target 0.95) ====\n")
cat("    [the slope is robust to IC: all methods ~nominal]\n")
show_cov(grepl("beta", agg$quantity))

# ---- figures --------------------------------------------------------------
cols <- c(oracle = "#1b9e77", MI = "#1f78b4", midpoint = "#e31a1c")
lab_task <- c(single = "Single event", cr = "Competing risks")
# facet_wrap (not facet_grid) so empty task x quantity combinations are dropped
fwrap <- function(...) facet_wrap(vars(task, quantity), ...,
  labeller = labeller(task = lab_task, quantity = label_parsed))

# Fig 1 (HEADLINE): coverage of FUNCTIONALS S(t)/CIF vs inspection-interval length
p_fun <- agg %>% filter(!grepl("beta", quantity)) %>%
  ggplot(aes(mean_int, coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  geom_line() + geom_point(size = 2) +
  fwrap(nrow = 1, scales = "free_x") +
  scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "mean inspection-interval length", y = "95% CI coverage (avg over t)",
    colour = NULL,
    title = "MI restores near-nominal coverage of S(t) and CIF; midpoint under-covers",
    subtitle = "degradation worsens as inspection intervals lengthen; MI tracks the oracle")
ggsave(file.path(FIGDIR, "fig1-functional-coverage.png"), p_fun,
  width = 10, height = 4.2, dpi = 130)

# Fig 2: coverage of the covariate TERMS (robust to IC under both methods)
p_term <- agg %>% filter(grepl("beta", quantity)) %>%
  ggplot(aes(mean_int, coverage, colour = method)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  geom_line() + geom_point(size = 2) +
  fwrap(nrow = 1) +
  scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "mean inspection-interval length", y = "95% CI coverage",
    colour = NULL,
    title = "Covariate-effect coverage is maintained by all methods",
    subtitle = "the regression slope is barely affected by interval censoring")
ggsave(file.path(FIGDIR, "fig2-term-coverage.png"), p_term,
  width = 10, height = 4.2, dpi = 130)

# Fig 3: point-estimate error (RMSE) of functionals vs inspection-interval length
p_rmse <- agg %>% filter(!grepl("beta", quantity)) %>%
  ggplot(aes(mean_int, rmse, colour = method)) +
  geom_line() + geom_point(size = 2) +
  fwrap(nrow = 1, scales = "free") +
  scale_colour_manual(values = cols) +
  labs(x = "mean inspection-interval length", y = "RMSE of point estimate",
    colour = NULL,
    title = "Point-estimate accuracy: MI <= midpoint, improving under heavy IC",
    subtitle = "both approach the oracle; gap to midpoint widens as intervals lengthen")
ggsave(file.path(FIGDIR, "fig3-point-rmse.png"), p_rmse,
  width = 10, height = 4.2, dpi = 130)

# Fig 4: functional CI width -- the MECHANISM behind the coverage gap:
# midpoint bands are too narrow (ignore imputation uncertainty); MI widens them.
p_w <- agg %>% filter(!grepl("beta", quantity)) %>%
  ggplot(aes(mean_int, width, colour = method)) +
  geom_line() + geom_point(size = 2) +
  fwrap(nrow = 1, scales = "free") +
  scale_colour_manual(values = cols) +
  labs(x = "mean inspection-interval length", y = "mean 95% CI width",
    colour = NULL,
    title = "Why midpoint under-covers: its S(t)/CIF bands are too narrow",
    subtitle = "MI widens the bands to reflect imputation uncertainty (oracle = exact times)")
ggsave(file.path(FIGDIR, "fig4-functional-width.png"), p_w,
  width = 10, height = 4.2, dpi = 130)

cat("\nFigures written to ", normalizePath(FIGDIR), "\n", sep = "")
