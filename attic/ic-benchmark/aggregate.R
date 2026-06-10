#!/usr/bin/env Rscript
# ===========================================================================
# IC benchmark: aggregation
# ---------------------------------------------------------------------------
# Reads results/raw/task_*.rds, defensively dedups, checks completeness, and
# writes pointwise + term-level summaries (with MC SEs) plus paired
# MI-vs-midpoint coverage differences to results/summary.rds / .csv.
# Usage: Rscript attic/ic-benchmark/aggregate.R
# ===========================================================================

suppressMessages(library(pammtools))
bench_dir <- if (file.exists("attic/ic-benchmark/config.R")) {
  "attic/ic-benchmark"
} else {
  "."
}
source(file.path(bench_dir, "config.R"))
source(file.path(bench_dir, "dgp.R"))
source(file.path(bench_dir, "metrics.R"))

files <- list.files(RAW_DIR, "^task_\\d+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no raw result files in ", RAW_DIR)
cat("reading", length(files), "raw files\n")
raw <- bind_rows(lapply(files, readRDS))
cat(nrow(raw), "rows\n")

# ---- defensive dedup ----------------------------------------------------------
# duplicate keys may arise from interrupted+resumed runs; deterministic columns
# must agree exactly (fit_time is non-deterministic and excluded)
key_cols <- c("cell_id", "rep", "method", "estimand", "t", "x")
det_cols <- c("est", "se", "lower", "upper", "truth", "ambiguity")
raw <- raw |>
  mutate(.key = do.call(paste, c(across(all_of(key_cols)), sep = "|")))
dup <- duplicated(raw$.key)
if (any(dup)) {
  conflict <- raw |>
    filter(.key %in% raw$.key[dup]) |>
    group_by(.key) |>
    summarise(n_val = n_distinct(across(all_of(det_cols))), .groups = "drop") |>
    filter(n_val > 1)
  if (nrow(conflict)) {
    stop(
      "dedup assertion failed: differing deterministic values for keys: ",
      paste(head(conflict$.key, 5), collapse = "; ")
    )
  }
  cat("dropping", sum(dup), "verified-identical duplicate rows\n")
  raw <- raw[!dup, ]
}
raw$.key <- NULL

# ---- completeness + failures -----------------------------------------------------
# compare against the configured rep count (env-aware), not the observed max,
# so a run truncated identically across all tasks is still flagged
expected <- TASK_TABLE |>
  filter(rep <= N_REP) |>
  distinct(cell_id, rep)
got <- raw |> distinct(cell_id, rep)
missing <- anti_join(expected, got, by = c("cell_id", "rep"))
if (nrow(missing)) {
  cat("WARNING:", nrow(missing), "missing (cell, rep) combinations:\n")
  print(count(missing, cell_id), n = 50)
}

failures <- raw |>
  group_by(cell_id, method) |>
  summarise(
    n_rep = n_distinct(rep),
    n_fail = n_distinct(rep[!is.na(error_msg)]),
    fail_rate = n_fail / n_rep,
    first_error = first(stats::na.omit(error_msg)),
    .groups = "drop"
  )
flagged <- filter(failures, fail_rate > 0.02)
if (nrow(flagged)) {
  cat("CELLS WITH >2% FAILURES (survivor-bias guard):\n")
  print(flagged)
}

# ---- pointwise summaries -----------------------------------------------------------
ok <- raw |> filter(is.na(error_msg), !is.na(estimand))

pointwise <- ok |>
  group_by(cell_id, method, estimand, t, x, misspecified) |>
  summarise(
    n_rep = sum(is.finite(err)),
    bias = mean(err, na.rm = TRUE),
    bias_mcse = stats::sd(err, na.rm = TRUE) / sqrt(n_rep),
    empSE = stats::sd(err, na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    avSE = sqrt(mean(se^2, na.rm = TRUE)),
    se_ratio = avSE / empSE,
    coverage = mean(covered, na.rm = TRUE),
    coverage_mcse = sqrt(coverage * (1 - coverage) / sum(!is.na(covered))),
    mean_width = mean(width, na.rm = TRUE),
    # bias-eliminated coverage: CI shifted to cover mean(est) instead of truth
    be_coverage = {
      mbar <- mean(est, na.rm = TRUE)
      mean(lower <= mbar & mbar <= upper, na.rm = TRUE)
    },
    # boundary-touching CIs only meaningful for probability/positive scales
    boundary_lo = ifelse(
      first(estimand) %in% c("surv", "hazard", "cif1", "cif2"),
      mean(lower <= 0, na.rm = TRUE),
      NA_real_
    ),
    boundary_hi = ifelse(
      first(estimand) %in% c("surv", "cif1", "cif2"),
      mean(upper >= 1, na.rm = TRUE),
      NA_real_
    ),
    mean_ambiguity = mean(ambiguity, na.rm = TRUE),
    .groups = "drop"
  )

# ---- term-level summaries (rep-clustered coverage MC SE) ----------------------------
term_level <- ok |>
  group_by(cell_id, method, estimand) |>
  group_modify(function(d, g) {
    cl <- mcse_coverage_clustered(d[!is.na(d$covered), ])
    rr <- avse_empse_ratio(d$se, d$err)
    tibble(
      n_rep = n_distinct(d$rep),
      bias = mean(d$err, na.rm = TRUE),
      bias_mcse = {
        per <- tapply(d$err, d$rep, mean, na.rm = TRUE)
        stats::sd(per, na.rm = TRUE) / sqrt(sum(is.finite(per)))
      },
      rmse = rmse_of(d$err),
      empSE = rr$empSE,
      avSE = rr$avSE,
      se_ratio = rr$ratio,
      coverage = cl$coverage,
      coverage_mcse = cl$mcse,
      mean_width = mean(d$width, na.rm = TRUE),
      mean_fit_time = mean(d$fit_time, na.rm = TRUE),
      misspecified = any(d$misspecified, na.rm = TRUE)
    )
  }) |>
  ungroup() |>
  left_join(
    CELLS |> select(cell_id, arm, baseline, effect, mechanism, rate, n, m),
    by = "cell_id"
  )

# per-cell censoring diagnostics: one value per rep (not per result row, which
# would weight reps by method row counts)
cens <- ok |>
  distinct(cell_id, rep, mean_int, prop_right_cens, prop_left_cens) |>
  group_by(cell_id) |>
  summarise(
    mean_int = mean(mean_int, na.rm = TRUE),
    prop_right_cens = mean(prop_right_cens, na.rm = TRUE),
    prop_left_cens = mean(prop_left_cens, na.rm = TRUE),
    .groups = "drop"
  )
term_level <- left_join(term_level, cens, by = "cell_id")

# ---- paired MI - midpoint coverage differences ----------------------------------------
paired <- ok |>
  filter(method %in% c("mi", "midpoint"), !is.na(covered)) |>
  group_by(cell_id, estimand) |>
  group_modify(function(d, g) {
    pd <- mcse_paired_coverage_diff(
      d[d$method == "mi", c("rep", "covered")],
      d[d$method == "midpoint", c("rep", "covered")]
    )
    tibble(diff = pd$diff, diff_mcse = pd$mcse, n_rep = pd$n_rep)
  }) |>
  ungroup() |>
  left_join(
    CELLS |> select(cell_id, arm, baseline, effect, rate, n),
    by = "cell_id"
  )

out <- list(
  pointwise = pointwise,
  term_level = term_level,
  paired = paired,
  failures = failures,
  missing = missing,
  cens = cens,
  n_raw_rows = nrow(raw)
)
saveRDS(out, file.path(RES_DIR, "summary.rds"))
write.csv(term_level, file.path(RES_DIR, "summary.csv"), row.names = FALSE)
cat("wrote", file.path(RES_DIR, "summary.rds"), "and summary.csv\n")
cat("\nterm-level preview:\n")
print(
  term_level |>
    filter(estimand %in% c("surv", "beta")) |>
    select(
      cell_id,
      method,
      estimand,
      coverage,
      coverage_mcse,
      rmse,
      se_ratio
    ) |>
    arrange(cell_id, estimand, method),
  n = 40
)
