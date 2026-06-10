#!/usr/bin/env Rscript
# ===========================================================================
# IC benchmark: local smoke / pilot driver
# ---------------------------------------------------------------------------
# Runs a subset of tasks locally by shelling into run-task.R semantics.
# Usage:
#   IC_SMOKE=1 Rscript attic/ic-benchmark/run-local.R           # smoke cells
#   IC_PILOT_REPS=20 Rscript attic/ic-benchmark/run-local.R 1 5 7   # tasks 1,5,7
# With no task arguments in smoke mode, one cheap cell per arm is selected.
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)

suppressMessages(library(pammtools))
bench_dir <- if (file.exists("attic/ic-benchmark/config.R")) {
  "attic/ic-benchmark"
} else {
  "."
}
source(file.path(bench_dir, "config.R"))

if (length(args) > 0) {
  task_ids <- as.integer(args)
} else if (IC_SMOKE) {
  # one cheap cell per arm (n smallest, rate densest within arm)
  task_ids <- CELLS |>
    group_by(arm) |>
    arrange(n, desc(rate), .by_group = TRUE) |>
    slice(1) |>
    pull(task_id)
} else {
  task_ids <- CELLS$task_id
}

n_cores <- max(1L, parallel::detectCores() - 1L)
cat(
  "run-local: tasks",
  paste(task_ids, collapse = ", "),
  "with",
  n_cores,
  "cores each\n"
)

rt <- file.path(bench_dir, "run-task.R")
for (id in task_ids) {
  status <- system2("Rscript", c(rt, id, n_cores))
  if (status != 0)
    cat(sprintf("run-local: task %d FAILED (exit %d)\n", id, status))
}
cat("run-local: done\n")
