#!/usr/bin/env Rscript
# ===========================================================================
# IC benchmark: run all reps of ONE cell (= one SLURM array task)
# ---------------------------------------------------------------------------
# Usage:  Rscript attic/ic-benchmark/run-task.R <task_id> [n_cores]
# Env:    IC_PILOT_REPS=20  -> pilot subset;  IC_SMOKE=1 -> tiny smoke settings
#
# Architecture: wave-batched mclapply. Forked children only RETURN result
# tibbles; the PARENT writes results/raw/task_<id>.rds after each wave via
# atomic tempfile + rename. Resume: completed rep ids in an existing file are
# skipped (their rows are kept).
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: run-task.R <task_id> [n_cores]")
task_id <- as.integer(args[1])
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1L

# isolated benchmark library (LRZ): must beat the .Rprofile-prepended shared
# lib, so prepend AFTER startup, right before loading
ic_lib <- Sys.getenv("IC_BENCH_LIB")
if (nzchar(ic_lib)) .libPaths(c(ic_lib, .libPaths()))
suppressMessages(library(pammtools))
stopifnot(exists("pamm_ic")) # guard against a shadowing pammtools install
# load namespaces in the PARENT so forked children don't each pay the cost
invisible(loadNamespace("icenReg"))
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

cell <- CELLS[CELLS$task_id == task_id, ]
stopifnot(nrow(cell) == 1)
reps <- TASK_TABLE |> filter(task_id == !!task_id, rep <= N_REP)
out_file <- file.path(RAW_DIR, sprintf("task_%02d.rds", task_id))

cat(sprintf(
  "[task %d] cell=%s n_rep=%d cores=%d nsim=%d%s\n",
  task_id,
  cell$cell_id,
  nrow(reps),
  n_cores,
  NSIM_RUN,
  if (IC_SMOKE) " (SMOKE)" else ""
))

# ---- resume -----------------------------------------------------------------
# rep-/fork-level failures are transient (OOM, fork death): retry them on
# resume; per-method errors from run_method() are deterministic and kept
done <- tibble()
if (file.exists(out_file)) {
  done <- readRDS(out_file)
  transient <- !is.na(done$error_msg) &
    grepl("^(REP|FORK)-LEVEL:", done$error_msg)
  retry_reps <- unique(done$rep[transient])
  done <- done |> filter(!rep %in% retry_reps)
  done_reps <- unique(done$rep)
  reps <- reps |> filter(!rep %in% done_reps)
  cat(sprintf(
    "[task %d] resume: %d reps done, %d transient retries, %d to go\n",
    task_id,
    length(done_reps),
    length(retry_reps),
    nrow(reps)
  ))
}

# ---- one replication ----------------------------------------------------------
run_rep <- function(rep_id, seed_data) {
  dat <- generate_data(cell, seed_data)
  subseeds <- method_subseeds(seed_data)
  methods <- get_applicable_methods(cell)
  res <- lapply(methods, function(meth) {
    run_method(meth, dat, cell, subseeds[METHOD_INDEX[meth]])
  })
  bind_rows(res) |>
    score_rep(cell) |>
    mutate(
      task_id = task_id,
      cell_id = cell$cell_id,
      rep = rep_id,
      seed_data = seed_data
    ) |>
    bind_cols(dat$meta)
}

# ---- wave-batched execution -----------------------------------------------------
atomic_save <- function(obj, file) {
  tmp <- paste0(file, ".tmp")
  saveRDS(obj, tmp)
  file.rename(tmp, file)
}

wave_size <- max(2L * n_cores, 8L)
waves <- split(seq_len(nrow(reps)), ceiling(seq_len(nrow(reps)) / wave_size))
t_start <- Sys.time()

for (w in seq_along(waves)) {
  ix <- waves[[w]]
  one <- function(i) {
    tryCatch(
      run_rep(reps$rep[i], reps$seed_data[i]),
      error = function(e) {
        # DGP-level failure (should be ~0); log with seed for replay
        tibble(
          estimand = NA_character_,
          t = NA_real_,
          x = NA_real_,
          est = NA_real_,
          se = NA_real_,
          lower = NA_real_,
          upper = NA_real_,
          ambiguity = NA_real_,
          error_msg = paste("REP-LEVEL:", conditionMessage(e)),
          method = NA_character_,
          fit_time = NA_real_,
          n_warnings = NA_integer_,
          first_warning = NA_character_,
          misspecified = NA,
          truth = NA_real_,
          covered = NA,
          width = NA_real_,
          err = NA_real_,
          task_id = task_id,
          cell_id = cell$cell_id,
          rep = reps$rep[i],
          seed_data = reps$seed_data[i],
          mean_int = NA_real_,
          prop_right_cens = NA_real_,
          prop_left_cens = NA_real_
        )
      }
    )
  }
  res <- if (n_cores > 1L) {
    parallel::mclapply(ix, one, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(ix, one)
  }
  # mclapply returns try-error objects on fork-level failure: convert, don't drop
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), TRUE)
  if (any(bad)) {
    for (j in which(bad)) {
      res[[j]] <- tibble(
        error_msg = paste("FORK-LEVEL:", as.character(res[[j]])),
        task_id = task_id,
        cell_id = cell$cell_id,
        rep = reps$rep[ix[j]],
        seed_data = reps$seed_data[ix[j]]
      )
    }
  }
  done <- bind_rows(done, bind_rows(res))
  atomic_save(done, out_file)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  n_done <- length(unlist(waves[seq_len(w)]))
  cat(sprintf(
    "[task %d] wave %d/%d done (%d/%d reps, %.1f s/rep)\n",
    task_id,
    w,
    length(waves),
    n_done,
    nrow(reps),
    el / n_done
  ))
}

cat(sprintf(
  "[task %d] complete: %d rows -> %s\n",
  task_id,
  nrow(done),
  out_file
))
