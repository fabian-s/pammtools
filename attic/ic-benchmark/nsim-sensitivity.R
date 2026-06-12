#!/usr/bin/env Rscript
# ===========================================================================
# IC benchmark: nsim sensitivity check (promised in PLAN P3)
# ---------------------------------------------------------------------------
# Coverage of the CI-bearing methods on ONE cell at nsim in {200, 500, 1000},
# same 20 reps/seeds. If coverage is stable across nsim, the production
# nsim = 500 is not a tail-order-statistic artifact.
# Writes results/nsim-sensitivity.rds (NOT results/raw -- different nsim rows
# would trip the dedup assertion there).
# Usage: Rscript attic/ic-benchmark/nsim-sensitivity.R [n_reps] [n_cores]
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 20L
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

cell <- CELLS[CELLS$cell_id == "core-peaked-ph-random-r0.6-n300-m10", ]
stopifnot(nrow(cell) == 1)
seeds <- TASK_TABLE |>
  filter(task_id == cell$task_id, rep <= n_reps)
methods <- c("mi", "midpoint", "oracle", "ic_par")

res <- list()
for (ns in c(200L, 500L, 1000L)) {
  NSIM_RUN <<- ns
  rows <- parallel::mclapply(
    seq_len(nrow(seeds)),
    function(i) {
      dat <- generate_data(cell, seeds$seed_data[i])
      sub <- method_subseeds(seeds$seed_data[i])
      bind_rows(lapply(methods, function(m) {
        run_method(m, dat, cell, sub[METHOD_INDEX[m]])
      })) |>
        score_rep(cell) |>
        mutate(rep = seeds$rep[i], nsim = ns)
    },
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  res[[as.character(ns)]] <- bind_rows(rows)
  cat("nsim", ns, "done\n")
}
res <- bind_rows(res)
saveRDS(res, file.path(RES_DIR, "nsim-sensitivity.rds"))

smry <- res |>
  filter(!is.na(covered)) |>
  group_by(nsim, method, estimand) |>
  group_modify(\(d, g) {
    cl <- mcse_coverage_clustered(d)
    tibble(coverage = cl$coverage, mcse = cl$mcse)
  }) |>
  ungroup() |>
  tidyr::pivot_wider(names_from = nsim, values_from = c(coverage, mcse))
print(as.data.frame(smry), digits = 3)
cat("wrote", file.path(RES_DIR, "nsim-sensitivity.rds"), "\n")
