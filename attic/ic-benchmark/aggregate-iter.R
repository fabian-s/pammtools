#!/usr/bin/env Rscript
# ===========================================================================
# Iterated-MI extension: aggregation
# ---------------------------------------------------------------------------
# Builds results/iter-summary.rds from results/iter-prod3 (mi_iter3, 24 core
# cells) + results/iter-prod5 (mi_iter5, 8 rate-0.3 cells) + the production
# mi rows (iter = 1), all paired rep-by-rep on identical data seeds.
# report.qmd reads iter-summary.rds only.
# Usage: Rscript attic/ic-benchmark/aggregate-iter.R
# ===========================================================================

ic_lib <- Sys.getenv("IC_BENCH_LIB")
if (nzchar(ic_lib)) .libPaths(c(ic_lib, .libPaths()))
suppressMessages(library(pammtools))
stopifnot(exists("pamm_ic"))
bench_dir <- if (file.exists("attic/ic-benchmark/config.R")) {
  "attic/ic-benchmark"
} else {
  "."
}
source(file.path(bench_dir, "config.R"))
source(file.path(bench_dir, "dgp.R"))
source(file.path(bench_dir, "metrics.R"))

read_dir <- function(sub) {
  fs <- list.files(
    file.path(RES_DIR, sub),
    "^task_\\d+\\.rds$",
    full.names = TRUE
  )
  if (!length(fs)) stop("no result files in ", sub)
  bind_rows(lapply(fs, readRDS))
}
it3 <- read_dir("iter-prod3")
it5 <- read_dir("iter-prod5")
cat(sprintf(
  "iter-prod3: %d cells, %d reps, %d failures | iter-prod5: %d cells, %d reps, %d failures\n",
  n_distinct(it3$cell_id),
  n_distinct(it3$rep),
  sum(!is.na(it3$error_msg)),
  n_distinct(it5$cell_id),
  n_distinct(it5$rep),
  sum(!is.na(it5$error_msg))
))

# matching production mi rows (the iter = 1 leg), same cells and reps
keys <- bind_rows(it3, it5) |> distinct(cell_id, rep)
prod_files <- list.files(RAW_DIR, "^task_\\d+\\.rds$", full.names = TRUE)
mi1 <- bind_rows(lapply(prod_files, function(f) {
  readRDS(f) |>
    filter(method == "mi", is.na(error_msg)) |>
    semi_join(keys, by = c("cell_id", "rep"))
}))

ok <- bind_rows(mi1, it3, it5) |>
  filter(is.na(error_msg), !is.na(estimand)) |>
  mutate(iter = dplyr::recode(method, mi = 1L, mi_iter3 = 3L, mi_iter5 = 5L))

# ---- term-level summaries (same construction as aggregate.R) -------------------
term_level <- ok |>
  group_by(cell_id, method, iter, estimand) |>
  group_modify(function(d, g) {
    cl <- mcse_coverage_clustered(d[!is.na(d$covered), ])
    tibble(
      n_rep = n_distinct(d$rep),
      bias = mean(d$err, na.rm = TRUE),
      bias_mcse = {
        per <- tapply(d$err, d$rep, mean, na.rm = TRUE)
        stats::sd(per, na.rm = TRUE) / sqrt(sum(is.finite(per)))
      },
      rmse = rmse_of(d$err),
      coverage = cl$coverage,
      coverage_mcse = cl$mcse,
      be_coverage = {
        bec <- d |>
          group_by(t, x) |>
          summarise(
            v = mean(
              lower <= mean(est, na.rm = TRUE) &
                mean(est, na.rm = TRUE) <= upper,
              na.rm = TRUE
            ),
            .groups = "drop"
          )
        mean(bec$v, na.rm = TRUE)
      },
      median_width = stats::median(d$width, na.rm = TRUE),
      mean_fit_time = mean(d$fit_time, na.rm = TRUE)
    )
  }) |>
  ungroup() |>
  left_join(
    CELLS |> select(cell_id, baseline, effect, rate, n),
    by = "cell_id"
  ) |>
  left_join(
    readRDS(file.path(RES_DIR, "summary.rds"))$cens |>
      select(cell_id, mean_int),
    by = "cell_id"
  )

# ---- pointwise (for the trajectory figure) --------------------------------------
pointwise <- ok |>
  filter(estimand %in% c("surv", "hazard")) |>
  group_by(cell_id, method, iter, estimand, t, x) |>
  summarise(
    truth = first(truth),
    bias = mean(err, na.rm = TRUE),
    bias_mcse = stats::sd(err, na.rm = TRUE) / sqrt(sum(is.finite(err))),
    coverage = mean(covered, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    CELLS |> select(cell_id, baseline, effect, rate, n),
    by = "cell_id"
  )

# ---- paired coverage differences vs one-step mi ----------------------------------
paired_iter <- ok |>
  filter(method != "mi", !is.na(covered)) |>
  group_by(method, iter, cell_id, estimand) |>
  group_modify(function(d, g) {
    base <- ok |>
      filter(
        method == "mi",
        cell_id == g$cell_id,
        estimand == g$estimand,
        !is.na(covered)
      )
    pd <- mcse_paired_coverage_diff(
      d[, c("rep", "covered")],
      base[, c("rep", "covered")]
    )
    tibble(diff = pd$diff, diff_mcse = pd$mcse, n_rep = pd$n_rep)
  }) |>
  ungroup() |>
  left_join(
    CELLS |> select(cell_id, baseline, effect, rate, n),
    by = "cell_id"
  )

out <- list(
  term_level = term_level,
  pointwise = pointwise,
  paired = paired_iter,
  failures = bind_rows(it3, it5) |>
    group_by(cell_id, method) |>
    summarise(
      n_rep = n_distinct(rep),
      n_fail = n_distinct(rep[!is.na(error_msg)]),
      .groups = "drop"
    )
)
saveRDS(out, file.path(RES_DIR, "iter-summary.rds"))
cat("wrote", file.path(RES_DIR, "iter-summary.rds"), "\n")
print(
  term_level |>
    filter(estimand == "surv") |>
    select(cell_id, method, coverage, coverage_mcse, bias) |>
    arrange(cell_id, method),
  n = 20
)
