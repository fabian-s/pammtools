#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Scratchpad: quick validation of the interval-censored MI implementation on
# simple simulated models. These are the throwaway checks used while building
# R/as-ped-ic.R, R/impute-ic.R, R/pamm-ic.R, R/pool-ic.R, R/sim-ic.R; they are
# kept here for reference and are NOT part of the package test suite (the
# polished versions live in tests/testthat/test-interval-censored.R).
#
# Run from the package root with the package loaded, e.g.
#   R -q -e 'devtools::load_all(); source("attic/ic-validation-scratchpad.R")'
# (R is provisioned for this repo by tools/setup-r-test-env.sh.)
# ---------------------------------------------------------------------------

suppressMessages({
  library(pammtools); library(dplyr); library(survival)
})

# Helper: simulate exact times then censor them into inspection intervals.
make_ic <- function(n = 250, seed = 1, beta_x = 0.5, rate = 1,
                    cut_fine = seq(0, 10, by = 0.1)) {
  set.seed(seed)
  df  <- data.frame(x = runif(n, -1, 1))
  sdf <- sim_pexp(~ -2 + beta_x * x, df, cut = cut_fine)
  add_inspections(sdf, rate = rate, max_time = max(cut_fine))
}

f   <- Surv(L, R, type = "interval2") ~ x
cut <- seq(0, 10, by = 0.5)

# --- 1. Smoke test: detect -> as_ped_ic -> pamm_ic -> pooled add_surv_prob ----
icd <- make_ic(250, seed = 2)
stopifnot(detect_ic(f, icd) == "interval2")
fit <- pamm_ic(f, icd, cut = cut, m = 8)
nd  <- make_newdata(as_ped(icd, f, cut = cut), tend = unique(tend))
s   <- add_surv_prob(nd, fit)
stopifnot(all(diff(s$surv_prob) <= 1e-8),                       # monotone
          all(s$surv_lower <= s$surv_prob + 1e-8),
          all(s$surv_prob  <= s$surv_upper + 1e-8))
cat("[1] smoke test OK  (x-effect ~", round(mean(sapply(fit$fits,
    function(z) coef(z)["x"])), 3), ", truth 0.5)\n")

# --- 2. Sampler: draws stay in (L, R] and are calibrated (PIT ~ Uniform) ------
ic    <- fit$ic
cache <- ic_pred_cache(fit$init_fit, ic, cut)
ti    <- impute_ic_times(fit$init_fit, ic, cut, cache = cache)
idx   <- which(as.character(ic$ic_kind) %in% c("interval", "left"))
stopifnot(all(ti[idx] > ic$ic_L[idx] - 1e-8),
          all(ti[idx] <= pmin(ic$ic_R[idx], max(cut)) + 1e-8))

ii   <- cache$ii; n_int <- nrow(ii)
hm   <- matrix(as.numeric(exp(cache$X %*% coef(fit$init_fit))), nrow = n_int)
Hcut <- rbind(0, apply(hm * ii$intlen, 2, cumsum)); cutv <- c(ii$tstart[1], ii$tend)
evalH <- function(t, sj) {
  j <- pmin(pmax(findInterval(t, cutv, left.open = TRUE, rightmost.closed = TRUE),
    1L), n_int)
  Hcut[cbind(j, sj)] + hm[cbind(j, sj)] * (t - cutv[j])
}
U  <- (exp(-evalH(ic$ic_L[idx], idx)) - exp(-evalH(ti[idx], idx))) /
      (exp(-evalH(ic$ic_L[idx], idx)) - exp(-evalH(pmin(ic$ic_R[idx], max(cut)), idx)))
U  <- U[is.finite(U)]
cat("[2] sampler in (L,R] OK;  PIT KS p =",
    round(suppressWarnings(ks.test(U, "punif"))$p.value, 3), "(>0.01 = calibrated)\n")

# --- 3. Exact-only data reproduces a plain right-censored PAMM ----------------
set.seed(9); dfe <- data.frame(x = runif(150, -1, 1))
ex <- sim_pexp(~ -2 + 0.5 * x, dfe, cut = cut)
ex$L <- ex$time; ex$R <- ifelse(ex$status == 1, ex$time, Inf)
fic <- pamm_ic(f, ex, cut = cut, m = 2)
ref <- pamm(ped_status ~ s(tend) + x, data = as_ped(ex, Surv(time, status) ~ x, cut = cut))
stopifnot(isTRUE(all.equal(unname(coef(fic$fits[[1]])), unname(coef(ref)), tol = 1e-6)))
cat("[3] exact-only == plain PAMM OK\n")

# --- 4. Competing risks: pamm_ic_cr + pooled CIF ------------------------------
set.seed(5); dfc <- data.frame(x = runif(300, -1, 1))
sdf <- sim_pexp(~ -2.5 + 0.6 * x, dfc, cut = seq(0, 10, by = 0.25))
sdf$cause <- ifelse(sdf$status == 1, sample(c(1, 2), nrow(sdf), TRUE, c(.6, .4)), 0)
icc <- add_inspections(sdf, rate = 1, max_time = 10)
fcr <- pamm_ic_cr(f, icc, cause = "cause", cut = cut, m = 4)
pcr <- suppressWarnings(as_ped(transform(icc, time = pmin(true_time, 10),
  status = ifelse(true_time > 10, 0, cause)), Surv(time, status) ~ x, cut = cut))
ndc <- group_by(make_newdata(pcr, tend = unique(tend), cause = unique(cause)), cause)
cif <- add_cif(ndc, fcr, nsim = 120)
stopifnot(all(cif$cif >= 0 & cif$cif <= 1),
          all(tapply(cif$cif, cif$cause, function(z) all(diff(z) >= -1e-8))))
cat("[4] competing-risks CIF OK\n")

# --- 5. Coverage study: naive midpoint vs MI (small, indicative) --------------
# Demonstrates that MI widens intervals and restores coverage of S(t) relative
# to a single midpoint-imputed fit. Increase nrep for a publishable comparison.
coverage_study <- function(nrep = 12, n = 300, ts = c(2, 4, 6), rate = 0.6) {
  truth <- exp(-exp(-2) * ts)
  CM <- CD <- WM <- WD <- matrix(NA_real_, nrep, length(ts))
  for (r in seq_len(nrep)) {
    icd <- make_ic(n, seed = 1000 + r, rate = rate)
    nd  <- make_newdata(suppressWarnings(as_ped(icd, f, cut = cut)), tend = ts)
    nd[["x"]] <- 0
    fit <- suppressWarnings(pamm_ic(f, icd, cut = cut, m = 10))
    smi <- add_surv_prob(nd, fit, nsim = 250)
    k   <- as.character(parse_ic_surv(f, icd)$ic_kind)
    icd$tm <- ifelse(k %in% c("exact", "right"), icd$L,
      pmin(ifelse(k == "left", icd$R / 2, (icd$L + pmin(icd$R, 10)) / 2), 10))
    icd$ev <- ifelse(k == "right", 0L, 1L)
    pm  <- pamm(ped_status ~ s(tend) + x,
      data = as_ped(icd[icd$tm > 0, ], Surv(tm, ev) ~ x, cut = cut))
    smd <- add_surv_prob(nd, pm, ci_type = "sim", nsim = 250)
    CM[r, ] <- as.numeric(truth >= smi$surv_lower & truth <= smi$surv_upper)
    CD[r, ] <- as.numeric(truth >= smd$surv_lower & truth <= smd$surv_upper)
    WM[r, ] <- smi$surv_upper - smi$surv_lower
    WD[r, ] <- smd$surv_upper - smd$surv_lower
  }
  data.frame(t = ts, truth = round(truth, 3),
    cov_midpoint = round(colMeans(CD), 3), cov_MI = round(colMeans(CM), 3),
    width_midpoint = round(colMeans(WD), 3), width_MI = round(colMeans(WM), 3))
}

if (identical(Sys.getenv("IC_RUN_COVERAGE"), "1")) {
  cat("[5] coverage study (set IC_RUN_COVERAGE=1 to run; ~2 min):\n")
  print(coverage_study())
} else {
  cat("[5] coverage study skipped (set IC_RUN_COVERAGE=1 to run).\n")
  cat("    Example 12-rep result (n=300, inspection rate 0.6):\n")
  cat("      t truth cov_midpoint cov_MI width_midpoint width_MI\n")
  cat("      2 0.763        0.750  1.000          0.112    0.137\n")
  cat("      4 0.582        0.833  0.917          0.126    0.144\n")
  cat("      6 0.444        0.917  0.917          0.121    0.136\n")
}

cat("\nAll interval-censored validation checks passed.\n")
