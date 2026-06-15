#!/usr/bin/env Rscript
# ===========================================================================
# Iterated-MI pilot analysis: mi_iter (iter = 2) vs production mi (one-step),
# paired rep-by-rep (identical data seeds). Run locally after fetching
# results/iter-pilot from LRZ.
# Usage: Rscript attic/ic-benchmark/compare-iter-pilot.R
# ===========================================================================

ic_lib <- Sys.getenv("IC_BENCH_LIB")
if (nzchar(ic_lib)) .libPaths(c(ic_lib, .libPaths()))
suppressMessages({
  library(dplyr)
  library(tidyr)
})
bench_dir <- if (file.exists("attic/ic-benchmark/config.R")) {
  "attic/ic-benchmark"
} else {
  "."
}

pilot_files <- list.files(
  file.path(bench_dir, "results", "iter-pilot"),
  "^task_\\d+\\.rds$",
  full.names = TRUE
)
if (!length(pilot_files)) stop("no iter-pilot result files")
it <- bind_rows(lapply(pilot_files, readRDS))
cat(sprintf(
  "iter-pilot: %d rows, %d cells, %d reps/cell, %d failures\n",
  nrow(it),
  n_distinct(it$cell_id),
  max(table(distinct(it, cell_id, rep)$cell_id)),
  sum(!is.na(it$error_msg))
))
it <- it |> filter(is.na(error_msg))

# matching production mi rows (same cells, same reps -> exactly paired data)
prod_files <- list.files(
  file.path(bench_dir, "results", "raw"),
  "^task_\\d+\\.rds$",
  full.names = TRUE
)
keys <- it |> distinct(cell_id, rep)
mi <- bind_rows(lapply(prod_files, function(f) {
  d <- readRDS(f)
  d |>
    filter(method == "mi", is.na(error_msg)) |>
    semi_join(keys, by = c("cell_id", "rep"))
}))
stopifnot(nrow(mi) > 0)

both <- bind_rows(it, mi) |>
  filter(estimand %in% c("surv", "hazard", "beta", "logHR"))

# ---- pointwise: the headline check is early-t survival bias -------------------
cat("\n== pointwise S(t|x=0): one-step vs iterated, paired ==\n")
ptw <- both |>
  filter(estimand == "surv", x == 0) |>
  select(cell_id, rep, method, t, err, covered) |>
  pivot_wider(names_from = method, values_from = c(err, covered))
ptw_sum <- ptw |>
  group_by(cell_id, t) |>
  summarise(
    n = sum(complete.cases(err_mi, err_mi_iter)),
    bias_mi = mean(err_mi, na.rm = TRUE),
    bias_iter = mean(err_mi_iter, na.rm = TRUE),
    # paired MC SE of the bias difference
    dbias = mean(err_mi_iter - err_mi, na.rm = TRUE),
    dbias_mcse = stats::sd(err_mi_iter - err_mi, na.rm = TRUE) / sqrt(n),
    cov_mi = mean(covered_mi, na.rm = TRUE),
    cov_iter = mean(covered_mi_iter, na.rm = TRUE),
    .groups = "drop"
  )
print(
  ptw_sum |>
    filter(t %in% c(1, 3, 7)) |>
    mutate(across(where(is.numeric), \(v) round(v, 3))) |>
    arrange(cell_id, t),
  n = 60
)

# ---- term-level coverage, paired -----------------------------------------------
cat("\n== term-level coverage (grid-avg), paired by rep ==\n")
tlc <- both |>
  filter(!is.na(covered)) |>
  group_by(cell_id, estimand, method, rep) |>
  summarise(c_rep = mean(covered), .groups = "drop") |>
  pivot_wider(names_from = method, values_from = c_rep) |>
  group_by(cell_id, estimand) |>
  summarise(
    n = sum(complete.cases(mi, mi_iter)),
    cov_mi = mean(mi, na.rm = TRUE),
    cov_iter = mean(mi_iter, na.rm = TRUE),
    diff = cov_iter - cov_mi,
    diff_mcse = stats::sd(mi_iter - mi, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )
print(
  tlc |> mutate(across(where(is.numeric), \(v) round(v, 3))),
  n = 60
)

# ---- timing ---------------------------------------------------------------------
cat("\n== fit time (s): mi_iter pilot vs production mi (same reps) ==\n")
tm <- both |>
  distinct(cell_id, rep, method, fit_time) |>
  group_by(cell_id, method) |>
  summarise(med = median(fit_time, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = method, values_from = med) |>
  mutate(ratio = mi_iter / mi)
print(tm |> mutate(across(where(is.numeric), \(v) round(v, 1))))

# ---- verdict ---------------------------------------------------------------------
cat("\n== verdict heuristics ==\n")
v1 <- ptw_sum |> filter(t == 1)
cat(sprintf(
  "S(1|0) bias: one-step %.3f-%.3f -> iterated %.3f-%.3f across cells; bias shrunk in %d/%d cells (significantly, |dbias| > 2 MC SE toward 0, in %d).\n",
  min(v1$bias_mi),
  max(v1$bias_mi),
  min(v1$bias_iter),
  max(v1$bias_iter),
  sum(abs(v1$bias_iter) < abs(v1$bias_mi)),
  nrow(v1),
  sum(abs(v1$bias_iter) < abs(v1$bias_mi) & abs(v1$dbias) > 2 * v1$dbias_mcse)
))
cat(sprintf(
  "Term-level surv coverage: improved in %d/%d cells (significantly in %d); largest gain %+.3f, largest loss %+.3f.\n",
  sum(tlc$diff[tlc$estimand == "surv"] > 0),
  sum(tlc$estimand == "surv"),
  sum(tlc$estimand == "surv" & tlc$diff > 2 * tlc$diff_mcse),
  max(tlc$diff[tlc$estimand == "surv"]),
  min(tlc$diff[tlc$estimand == "surv"])
))
