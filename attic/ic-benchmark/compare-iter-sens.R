#!/usr/bin/env Rscript
# ===========================================================================
# iter-sensitivity analysis: bias/coverage trajectory over iter in {1, 2, 3, 5}
# (iter 1 = production mi, iter 2 = iter-pilot, iter 3/5 = iter-sens), all
# paired rep-by-rep on identical data seeds.
# Usage: Rscript attic/ic-benchmark/compare-iter-sens.R
# ===========================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
})
bench_dir <- if (file.exists("attic/ic-benchmark/config.R")) {
  "attic/ic-benchmark"
} else {
  "."
}

read_dir <- function(sub) {
  fs <- list.files(
    file.path(bench_dir, "results", sub),
    "^task_\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(fs)) stop("no result files in ", sub)
  bind_rows(lapply(fs, readRDS))
}
sens <- read_dir("iter-sens") |> filter(is.na(error_msg))
it2 <- read_dir("iter-pilot") |> filter(is.na(error_msg))
keys <- sens |> distinct(cell_id, rep)
it2 <- semi_join(it2, keys, by = c("cell_id", "rep"))
prod_files <- list.files(
  file.path(bench_dir, "results", "raw"),
  "^task_\\d+\\.rds$",
  full.names = TRUE
)
mi1 <- bind_rows(lapply(prod_files, function(f) {
  readRDS(f) |>
    filter(method == "mi", is.na(error_msg)) |>
    semi_join(keys, by = c("cell_id", "rep"))
}))

all_it <- bind_rows(mi1, it2, sens) |>
  mutate(
    iter = recode(method, mi = 1L, mi_iter = 2L, mi_iter3 = 3L, mi_iter5 = 5L)
  ) |>
  filter(estimand %in% c("surv", "hazard"))
cat(sprintf(
  "rows: %d | cells: %d | reps/cell: %d | iters: %s | failures already dropped\n",
  nrow(all_it),
  n_distinct(all_it$cell_id),
  max(table(distinct(all_it, cell_id, rep)$cell_id)),
  paste(sort(unique(all_it$iter)), collapse = ",")
))

# ---- trajectory: S(t|0) bias and coverage by iter ------------------------------
cat("\n== S(t|0) bias / coverage trajectory over iter ==\n")
tr <- all_it |>
  filter(estimand == "surv", x == 0, t %in% c(1, 2, 3)) |>
  group_by(cell_id, t, iter) |>
  summarise(
    bias = mean(err, na.rm = TRUE),
    bias_mcse = stats::sd(err, na.rm = TRUE) / sqrt(sum(is.finite(err))),
    coverage = mean(covered, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(cell_id, t, iter)
print(
  tr |>
    pivot_wider(
      id_cols = c(cell_id, t),
      names_from = iter,
      values_from = c(bias, coverage)
    ) |>
    mutate(across(where(is.numeric), \(v) round(v, 3))),
  n = 30,
  width = Inf
)

# ---- term-level coverage trajectory ----------------------------------------------
cat("\n== term-level coverage trajectory (surv / hazard) ==\n")
tl_tr <- all_it |>
  filter(!is.na(covered)) |>
  group_by(cell_id, estimand, iter) |>
  summarise(coverage = mean(covered, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = iter, values_from = coverage, names_prefix = "it") |>
  mutate(across(where(is.numeric), \(v) round(v, 3)))
print(tl_tr, n = 20)

# ---- fit-time scaling -------------------------------------------------------------
cat("\n== median fit time (s) by iter ==\n")
print(
  all_it |>
    distinct(cell_id, rep, iter, fit_time) |>
    group_by(cell_id, iter) |>
    summarise(
      med = round(median(fit_time, na.rm = TRUE), 1),
      .groups = "drop"
    ) |>
    pivot_wider(names_from = iter, values_from = med, names_prefix = "it")
)

# ---- verdict ---------------------------------------------------------------------
cat("\n== verdict ==\n")
v <- tr |> filter(t == 1)
for (cl in unique(v$cell_id)) {
  d <- v[v$cell_id == cl, ]
  cat(sprintf(
    "%s: S(1|0) bias %s (MC SE ~%.3f); coverage %s\n",
    cl,
    paste(sprintf("it%d %+.3f", d$iter, d$bias), collapse = ", "),
    max(d$bias_mcse),
    paste(sprintf("%.2f", d$coverage), collapse = " -> ")
  ))
}
conv <- v |>
  group_by(cell_id) |>
  summarise(
    floor_hit = abs(bias[iter == 5]) > 2 * bias_mcse[iter == 5],
    gain_3_5 = abs(bias[iter == 3]) - abs(bias[iter == 5]),
    .groups = "drop"
  )
cat(sprintf(
  "\nBias floor (S(1|0) bias at iter 5 still > 2 MC SE from 0): %d/%d cells; mean |bias| reduction from iter 3 -> 5: %.3f.\n",
  sum(conv$floor_hit),
  nrow(conv),
  mean(conv$gain_3_5)
))
