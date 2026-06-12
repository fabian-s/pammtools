# ===========================================================================
# IC benchmark: scoring + Monte Carlo SE helpers
# ---------------------------------------------------------------------------
# score_rep(): join method estimates to the cell truth -> per-row metrics.
# Aggregation (with rep-clustered + paired MC SEs) lives in aggregate.R.
# Requires config.R + dgp.R sourced first.
# ===========================================================================

# Join one rep's stacked method estimates to the cell truth and compute
# per-row covered/width/err. `est_df` has columns from run_method():
#   estimand, t, x, est, se, lower, upper, ambiguity, method, fit_time,
#   n_warnings, first_warning, misspecified, error_msg
score_rep <- function(est_df, cell) {
  truth <- make_truth_table(cell) # estimand, t, x, truth
  # robust key join (t/x may be NA for beta/logHR)
  # as.character is elementwise (no vector-dependent padding like format())
  key <- function(d) {
    paste(d$estimand, as.character(d$t), as.character(d$x), sep = "|")
  }
  est_df$.k <- key(est_df)
  truth$.k <- key(truth)
  stopifnot(!anyDuplicated(truth$.k))
  out <- est_df
  out$truth <- truth$truth[match(out$.k, truth$.k)]
  # every non-error estimate row must match a truth row -- a silent mismatch
  # (typo, float key, unexpected estimand) would vanish from all summaries
  unmatched <- is.na(out$truth) & !is.na(out$estimand) & is.na(out$error_msg)
  if (any(unmatched)) {
    stop(
      "score_rep: ",
      sum(unmatched),
      " estimate rows without truth match, e.g. ",
      paste(head(out$.k[unmatched], 3), collapse = "; ")
    )
  }
  out$.k <- NULL
  out$covered <- out$lower <= out$truth & out$truth <= out$upper
  out$width <- out$upper - out$lower
  out$err <- out$est - out$truth
  out
}

# ---- Monte Carlo SE helpers (used in aggregate.R) ----------------------------
# Term-level coverage clusters indicators at the rep level: average over the
# (t, x) grid WITHIN each rep, then take sd across reps / sqrt(n_rep).
# d: rows for one (cell, method, estimand) group with columns rep, covered.
mcse_coverage_clustered <- function(d) {
  per_rep <- tapply(d$covered, d$rep, mean, na.rm = TRUE)
  per_rep <- per_rep[is.finite(per_rep)]
  n <- length(per_rep)
  if (n == 0) {
    return(list(coverage = NA_real_, mcse = NA_real_, n_rep = 0L))
  }
  list(coverage = mean(per_rep), mcse = stats::sd(per_rep) / sqrt(n), n_rep = n)
}

# Pointwise coverage MC SE (single t,x): binomial.
mcse_coverage_binom <- function(covered) {
  covered <- covered[!is.na(covered)]
  n <- length(covered)
  p <- mean(covered)
  list(coverage = p, mcse = sqrt(p * (1 - p) / n), n = n)
}

# Bias + its MC SE = empSE / sqrt(n).
mcse_bias <- function(err) {
  err <- err[is.finite(err)]
  n <- length(err)
  list(
    bias = mean(err),
    mcse = stats::sd(err) / sqrt(n),
    empSE = stats::sd(err),
    n = n
  )
}

# RMSE.
rmse_of <- function(err) sqrt(mean(err^2, na.rm = TRUE))

# Paired coverage difference (method A vs B on the SAME reps): rep-clustered.
# dA, dB: data frames with rep, covered for one estimand; aligned on rep.
mcse_paired_coverage_diff <- function(dA, dB) {
  a <- tapply(dA$covered, dA$rep, mean, na.rm = TRUE)
  b <- tapply(dB$covered, dB$rep, mean, na.rm = TRUE)
  common <- intersect(names(a), names(b))
  diff <- a[common] - b[common]
  diff <- diff[is.finite(diff)]
  n <- length(diff)
  list(diff = mean(diff), mcse = stats::sd(diff) / sqrt(n), n_rep = n)
}

# avSE/empSE ratio (model-based SE vs empirical SE): >1 over-, <1 under-states.
avse_empse_ratio <- function(se, err) {
  ok <- is.finite(se) & is.finite(err)
  list(
    avSE = sqrt(mean(se[ok]^2)),
    empSE = stats::sd(err[ok]),
    ratio = sqrt(mean(se[ok]^2)) / stats::sd(err[ok])
  )
}
