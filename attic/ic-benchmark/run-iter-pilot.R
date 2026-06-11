#!/usr/bin/env Rscript
# ===========================================================================
# IC benchmark EXTENSION pilot: iterated MI (pamm_ic iter = 2)
# ---------------------------------------------------------------------------
# Runs ONLY the mi_iter method on the cells where the one-step early-t bias
# was found, using the SAME per-rep data seeds as production, so every rep
# pairs exactly with the existing mi/midpoint/oracle production results.
# Writes results/iter-pilot/task_<k>.rds -- production results are untouched.
#
# Usage:  Rscript attic/ic-benchmark/run-iter-pilot.R <pilot_task> [n_cores]
#         pilot_task indexes ITER_PILOT_CELLS (1-10)
# Env:    ITER_PILOT_REPS    (default 50)
#         ITER_PILOT_METHODS (default "mi_iter"; comma-separated, e.g.
#                             "mi_iter3,mi_iter5" for the iter-sensitivity run)
#         ITER_PILOT_OUTDIR  (default "iter-pilot"; results/<outdir>/)
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: run-iter-pilot.R <pilot_task> [n_cores]")
pilot_task <- as.integer(args[1])
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1L

ic_lib <- Sys.getenv("IC_BENCH_LIB")
if (nzchar(ic_lib)) .libPaths(c(ic_lib, .libPaths()))
suppressMessages(library(pammtools))
stopifnot(exists("pamm_ic")) # guard against a shadowing pammtools install
stopifnot("iter" %in% names(formals(pamm_ic))) # extension must be installed
invisible(loadNamespace("mvtnorm"))
bench_dir <- if (file.exists("attic/ic-benchmark/config.R")) {
  "attic/ic-benchmark"
} else {
  "."
}
source(file.path(bench_dir, "config.R"))
source(file.path(bench_dir, "dgp.R"))
source(file.path(bench_dir, "methods.R"))
source(file.path(bench_dir, "metrics.R"))

# the 8 sparsest core cells (one-step bias largest) + the two moderate-rate
# large-n cells where S coverage already dipped (0.819/0.821)
ITER_PILOT_CELLS <- c(
  "core-const-ph-random-r0.3-n300-m10",
  "core-const-ph-random-r0.3-n1000-m10",
  "core-peaked-ph-random-r0.3-n300-m10",
  "core-peaked-ph-random-r0.3-n1000-m10",
  "core-const-tv-random-r0.3-n300-m10",
  "core-const-tv-random-r0.3-n1000-m10",
  "core-peaked-tv-random-r0.3-n300-m10",
  "core-peaked-tv-random-r0.3-n1000-m10",
  "core-const-ph-random-r0.6-n1000-m10",
  "core-peaked-ph-random-r0.6-n1000-m10"
)
stopifnot(all(ITER_PILOT_CELLS %in% CELLS$cell_id))
stopifnot(pilot_task >= 1, pilot_task <= length(ITER_PILOT_CELLS))

n_pilot_reps <- {
  v <- Sys.getenv("ITER_PILOT_REPS")
  if (nzchar(v)) as.integer(v) else 50L
}
methods_run <- strsplit(
  Sys.getenv("ITER_PILOT_METHODS", "mi_iter"),
  ","
)[[1]]
stopifnot(all(methods_run %in% names(METHOD_INDEX)))
out_sub <- Sys.getenv("ITER_PILOT_OUTDIR", "iter-pilot")

cell <- CELLS[CELLS$cell_id == ITER_PILOT_CELLS[pilot_task], ]
stopifnot(nrow(cell) == 1)
# SAME seeds as production: reps 1..n_pilot_reps of this cell's task
reps <- TASK_TABLE |>
  filter(task_id == cell$task_id, rep <= n_pilot_reps)

out_dir <- file.path(RES_DIR, out_sub)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("task_%02d.rds", pilot_task))

cat(sprintf(
  "[iter-pilot %d] cell=%s methods=%s n_rep=%d cores=%d\n",
  pilot_task,
  cell$cell_id,
  paste(methods_run, collapse = "+"),
  nrow(reps),
  n_cores
))

run_rep <- function(rep_id, seed_data) {
  dat <- generate_data(cell, seed_data)
  sub <- method_subseeds(seed_data)
  bind_rows(lapply(methods_run, function(meth) {
    run_method(meth, dat, cell, sub[METHOD_INDEX[meth]]) |>
      score_rep(cell)
  })) |>
    mutate(
      task_id = cell$task_id,
      cell_id = cell$cell_id,
      rep = rep_id,
      seed_data = seed_data
    ) |>
    bind_cols(dat$meta)
}

# resume: keep completed reps
done <- tibble()
if (file.exists(out_file)) {
  done <- readRDS(out_file)
  reps <- reps |> filter(!rep %in% unique(done$rep))
  cat(sprintf("resume: %d reps already done\n", length(unique(done$rep))))
}

atomic_save <- function(obj, file) {
  tmp <- paste0(file, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, file)
}

wave_size <- max(2L * n_cores, 8L)
waves <- split(seq_len(nrow(reps)), ceiling(seq_len(nrow(reps)) / wave_size))
t0 <- Sys.time()
for (w in seq_along(waves)) {
  ix <- waves[[w]]
  one <- function(i) {
    tryCatch(
      run_rep(reps$rep[i], reps$seed_data[i]),
      error = function(e) {
        tibble(
          rep = reps$rep[i],
          seed_data = reps$seed_data[i],
          cell_id = cell$cell_id,
          method = "mi_iter",
          error_msg = paste("REP-LEVEL:", conditionMessage(e))
        )
      }
    )
  }
  res <- if (n_cores > 1L) {
    parallel::mclapply(ix, one, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(ix, one)
  }
  bad <- vapply(res, inherits, TRUE, what = "try-error")
  if (any(bad)) {
    res[bad] <- lapply(which(bad), function(i) {
      tibble(
        rep = reps$rep[ix[i]],
        cell_id = cell$cell_id,
        method = "mi_iter",
        error_msg = "FORK-LEVEL: mclapply try-error"
      )
    })
  }
  done <- bind_rows(done, bind_rows(res))
  atomic_save(done, out_file)
  cat(sprintf(
    "wave %d/%d done (%d reps, %.1f min elapsed)\n",
    w,
    length(waves),
    length(unique(done$rep)),
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}
cat(sprintf(
  "[iter-pilot %d] DONE: %d reps, %d failures\n",
  pilot_task,
  length(unique(done$rep)),
  sum(!is.na(done$error_msg))
))
