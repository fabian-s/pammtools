#!/usr/bin/env Rscript
# ===========================================================================
# IC benchmark: Phase 0 known-answer validation
# ---------------------------------------------------------------------------
# Validates the INFRASTRUCTURE on problems with known answers before any
# benchmark code is trusted:
#  (a) exact truth vs Kaplan-Meier of 1e6 exact draws (max |diff| < 0.005)
#  (b) oracle PAMM on const+PH, n=1000, 200 reps: beta & S(t) coverage within
#      +-3 MC SE of 0.95 (diagnostic gate)
#  (c) pamm_ic with near-exact intervals (rate=20) ~ oracle on the same reps
#  (d) Turnbull cross-check: icenReg::ic_np vs survival::survfit NPMLE
#  (e) icenReg API checks (covarOffset, getFitEsts, getSCurves shapes)
# Usage: Rscript attic/ic-benchmark/phase0-known-answer.R [n_reps] [n_cores]
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 200L
n_cores <- if (length(args) >= 2) {
  as.integer(args[2])
} else {
  max(1L, parallel::detectCores() - 1L)
}

ic_lib <- Sys.getenv("IC_BENCH_LIB")
if (nzchar(ic_lib)) .libPaths(c(ic_lib, .libPaths()))
suppressMessages(library(pammtools))
stopifnot(exists("pamm_ic")) # guard against a shadowing pammtools install
bench_dir <- if (file.exists("attic/ic-benchmark/config.R")) {
  "attic/ic-benchmark"
} else {
  "."
}
source(file.path(bench_dir, "config.R"))
source(file.path(bench_dir, "dgp.R"))
source(file.path(bench_dir, "methods.R"))
source(file.path(bench_dir, "metrics.R"))

ok_all <- TRUE
check <- function(cond, msg) {
  status <- if (isTRUE(cond)) "ok " else "FAIL"
  if (!isTRUE(cond)) ok_all <<- FALSE
  cat(sprintf("[%s] %s\n", status, msg))
}

cell_p0 <- CELLS[CELLS$cell_id == "core-const-ph-random-r1.5-n1000-m10", ]
stopifnot(nrow(cell_p0) == 1)

# ---- (a) exact truth vs KM of many exact draws ---------------------------------
# sim_pexp splits internally on FINE (200 intervals/subject): simulate in
# chunks to bound memory. 4e5 draws/stratum -> KM MC SE ~ 0.0016 at S = 0.5.
cat("\n--- (a) truth vs KM (chunked large-sample draws) ---\n")
sim_times_chunked <- function(form, x_value, n_total, chunk = 5e4, seed) {
  set.seed(seed)
  out <- vector("list", ceiling(n_total / chunk))
  for (k in seq_along(out)) {
    nk <- min(chunk, n_total - (k - 1) * chunk)
    d <- tibble(id = seq_len(nk), x_num = rep(x_value, nk)) |>
      mutate(x = factor(x_num, levels = c(0, 1)))
    s <- sim_pexp(form, d, cut = FINE)
    out[[k]] <- s[, c("time", "status")]
    rm(s, d)
  }
  bind_rows(out)
}
tt <- make_truth_table(cell_p0)
for (xv in c(0, 1)) {
  sdf <- sim_times_chunked(cell_sim_formula(cell_p0), xv, 4e5, seed = 101 + xv)
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = sdf)
  km_at <- summary(km, times = TEVAL)$surv
  tr <- tt$truth[tt$estimand == "surv" & tt$x == xv]
  dmax <- max(abs(km_at - tr))
  check(
    dmax < 0.005,
    sprintf("const x=%d: max |KM - truth| = %.4f < 0.005", xv, dmax)
  )
  rm(sdf)
  invisible(gc())
}

# also peaked-baseline truth (the harder case)
cell_pk <- CELLS[CELLS$cell_id == "core-peaked-ph-random-r1.5-n1000-m10", ]
sdf <- sim_times_chunked(cell_sim_formula(cell_pk), 0L, 4e5, seed = 103)
km <- survival::survfit(survival::Surv(time, status) ~ 1, data = sdf)
km_at <- summary(km, times = TEVAL)$surv
tr <- make_truth_table(cell_pk)
tr <- tr$truth[tr$estimand == "surv" & tr$x == 0]
dmax <- max(abs(km_at - tr))
check(
  dmax < 0.005,
  sprintf("peaked x=0: max |KM - truth| = %.4f < 0.005", dmax)
)
rm(sdf)
invisible(gc())

# ---- (b)+(c) coverage of oracle and near-exact pamm_ic -------------------------
cat(sprintf(
  "\n--- (b)+(c) oracle & near-exact MI coverage (%d reps) ---\n",
  n_reps
))
cell_ne <- cell_p0
cell_ne$rate <- 20 # near-exact inspection intervals

one_rep <- function(i) {
  seed <- 7e6 + i
  dat <- generate_data(cell_p0, seed)
  sub <- method_subseeds(seed)
  ora <- run_method("oracle", dat, cell_p0, sub[METHOD_INDEX["oracle"]]) |>
    score_rep(cell_p0)
  dat_ne <- generate_data(cell_ne, seed)
  mi <- run_method("mi", dat_ne, cell_ne, sub[METHOD_INDEX["mi"]]) |>
    score_rep(cell_ne)
  bind_rows(
    mutate(ora, method = "oracle"),
    mutate(mi, method = "mi_near_exact")
  ) |>
    mutate(rep = i)
}
res <- parallel::mclapply(
  seq_len(n_reps),
  function(i) {
    tryCatch(one_rep(i), error = function(e) {
      tibble(rep = i, error_msg = conditionMessage(e))
    })
  },
  mc.cores = n_cores,
  mc.preschedule = FALSE
)
res <- bind_rows(res)
errs <- res |> filter(!is.na(error_msg)) |> distinct(rep, error_msg)
if (nrow(errs)) {
  cat("rep-level errors:\n")
  print(head(errs, 10))
}
n_eff <- length(setdiff(unique(res$rep), errs$rep))
check(
  n_eff >= 0.98 * n_reps,
  sprintf("%d/%d reps completed (>= 98%% required)", n_eff, n_reps)
)

gate <- function(d, label, target = 0.95) {
  cl <- mcse_coverage_clustered(d)
  lo <- target - 3 * sqrt(target * (1 - target) / cl$n_rep)
  hi <- target + 3 * sqrt(target * (1 - target) / cl$n_rep)
  check(
    cl$coverage >= lo & cl$coverage <= hi,
    sprintf(
      "%s coverage %.3f (MC SE %.3f) in [%.3f, %.3f]",
      label,
      cl$coverage,
      cl$mcse,
      lo,
      hi
    )
  )
}
for (meth in c("oracle", "mi_near_exact")) {
  d <- res[res$method == meth & !is.na(res$covered), ]
  gate(d[d$estimand == "surv", ], paste(meth, "S(t|x)"))
  gate(d[d$estimand == "beta", ], paste(meth, "beta"))
  gate(d[d$estimand == "hazard", ], paste(meth, "hazard"))
}

# ---- (d) Turnbull cross-check: ic_np vs survival::survfit ----------------------
cat("\n--- (d) Turnbull: ic_np vs survfit ---\n")
set.seed(103)
dat <- generate_data(
  CELLS[CELLS$cell_id == "core-peaked-ph-random-r0.6-n300-m10", ],
  555L
)
icd0 <- dat$icd[dat$icd$x_num == 0, ]
np <- icenReg::ic_np(cbind(L, R) ~ 0, data = icd0[, c("L", "R")])
s_np <- turnbull_S_at(np, TEVAL)
sf <- survival::survfit(
  survival::Surv(L, R, type = "interval2") ~ 1,
  data = icd0
)
# survfit NPMLE: evaluate its step function at TEVAL
sf_fun <- stats::stepfun(sf$time, c(1, sf$surv))
s_sf <- sf_fun(TEVAL)
cat("ic_np (midpoint convention):", round(s_np$est, 3), "\n")
cat("survfit:                    ", round(s_sf, 3), "\n")
# agreement up to the Turnbull-interval ambiguity
check(
  all(abs(s_np$est - s_sf) <= s_np$ambiguity + 0.02),
  "ic_np and survfit agree within Turnbull ambiguity + 0.02"
)

# ---- (e) icenReg API invariants -------------------------------------------------
cat("\n--- (e) icenReg API invariants ---\n")
fit <- icenReg::ic_par(
  survival::Surv(L, R, type = "interval2") ~ x_num,
  data = dat$icd,
  model = "ph",
  dist = "weibull"
)
check(
  all(c("log_shape", "log_scale", "x_num") %in% names(fit$coefficients)),
  "ic_par coefficient names"
)
check(is.numeric(as.numeric(fit$covarOffset)), "ic_par covarOffset present")
off <- as.numeric(fit$covarOffset)
cf <- fit$coefficients
S_man <- exp(-exp(exp(cf[1]) * (log(TEVAL) - cf[2])) * exp(cf[3] * (1 - off)))
S_pkg <- 1 -
  icenReg::getFitEsts(fit, newdata = data.frame(x_num = 1), q = TEVAL)
check(max(abs(S_man - S_pkg)) < 1e-6, "ic_par closed form == getFitEsts (x=1)")

sp <- icenReg::ic_sp(
  survival::Surv(L, R, type = "interval2") ~ x_num,
  data = dat$icd,
  bs_samples = 10
)
check(
  is.matrix(sp$bsMat) && "x_num" %in% colnames(sp$bsMat),
  "ic_sp bsMat has coefficient column"
)

cat(sprintf(
  "\n%s\n",
  if (ok_all) "PHASE 0: ALL CHECKS PASSED" else "PHASE 0: FAILURES"
))
quit(status = as.integer(!ok_all))
