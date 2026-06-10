#!/usr/bin/env Rscript
# ===========================================================================
# IC benchmark: P1 single-fit debugging
# ---------------------------------------------------------------------------
# Runs every applicable method on ONE rep of selected cells, verbose, and
# checks schema conformance + method-specific invariants:
#   - all extractors return the standardized columns
#   - estimates join to truth (no NA truth rows except error rows)
#   - S(t|x) monotone non-increasing within stratum (grouped-grid guard)
#   - PH cell: logHR(t) flat == beta-hat for PH-model methods
#   - ic_par closed forms cross-checked against icenReg::getFitEsts
# Usage: Rscript attic/ic-benchmark/debug-single-fit.R [task_id] [rep]
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)

suppressMessages(library(pammtools))
bench_dir <- if (file.exists("attic/ic-benchmark/config.R")) {
  "attic/ic-benchmark"
} else {
  "."
}
source(file.path(bench_dir, "config.R"))
source(file.path(bench_dir, "dgp.R"))
source(file.path(bench_dir, "methods.R"))
source(file.path(bench_dir, "metrics.R"))

SCHEMA_COLS <- c(
  "estimand",
  "t",
  "x",
  "est",
  "se",
  "lower",
  "upper",
  "ambiguity",
  "error_msg",
  "method",
  "fit_time",
  "n_warnings",
  "first_warning",
  "misspecified",
  "truth",
  "covered",
  "width",
  "err"
)

debug_cells <- if (length(args) >= 1) {
  as.integer(args[1])
} else {
  # default: one PH core cell (all 6 methods), one TV cell, the icsp-boot
  # variant, one smooth and one CR-unknown cell
  c(
    CELLS$task_id[CELLS$cell_id == "core-peaked-ph-random-r0.6-n300-m10"],
    CELLS$task_id[CELLS$cell_id == "core-peaked-tv-random-r0.6-n300-m10"],
    CELLS$task_id[CELLS$cell_id == "core-peaked-ph-random-r1.5-n300-m10"],
    CELLS$task_id[CELLS$arm == "smooth"][1],
    CELLS$task_id[CELLS$arm == "cr_unknown"]
  )
}
rep_id <- if (length(args) >= 2) as.integer(args[2]) else 1L

ok_all <- TRUE
check <- function(cond, msg) {
  status <- if (isTRUE(cond)) "ok " else "FAIL"
  if (!isTRUE(cond)) ok_all <<- FALSE
  cat(sprintf("  [%s] %s\n", status, msg))
}

for (tid in debug_cells) {
  cell <- CELLS[CELLS$task_id == tid, ]
  seed <- TASK_TABLE$seed_data[
    TASK_TABLE$task_id == tid & TASK_TABLE$rep == rep_id
  ]
  cat(sprintf(
    "\n== cell %s (task %d, rep %d, seed %d) ==\n",
    cell$cell_id,
    tid,
    rep_id,
    seed
  ))
  dat <- generate_data(cell, seed)
  subseeds <- method_subseeds(seed)
  for (meth in get_applicable_methods(cell)) {
    cat(sprintf("-- %s --\n", meth))
    res <- run_method(meth, dat, cell, subseeds[METHOD_INDEX[meth]]) |>
      score_rep(cell)
    check(all(SCHEMA_COLS %in% names(res)), "schema columns present")
    check(all(is.na(res$error_msg)), paste("no errors:", res$error_msg[1]))
    if (all(is.na(res$error_msg))) {
      check(!any(is.na(res$truth)), "all rows joined to truth")
      sv <- res[res$estimand == "surv", ]
      if (nrow(sv) > 0) {
        mono <- all(tapply(sv$est, sv$x, function(s) all(diff(s) <= 1e-10)))
        check(mono, "S(t|x) monotone within stratum")
        check(all(sv$est >= 0 & sv$est <= 1), "S in [0,1]")
      }
      lh <- res[res$estimand == "logHR", ]
      if (
        nrow(lh) > 0 &&
          cell$model_type == "ph" &&
          meth %in% c("mi", "midpoint", "oracle", "ic_par", "ic_sp")
      ) {
        check(diff(range(lh$est)) < 1e-8, "PH model: logHR(t) flat")
        bt <- res[res$estimand == "beta", ]
        if (nrow(bt) == 1) {
          check(abs(lh$est[1] - bt$est) < 1e-6, "logHR == beta-hat (PH)")
        }
      }
      print(
        as.data.frame(
          res[, c(
            "estimand",
            "t",
            "x",
            "est",
            "lower",
            "upper",
            "truth",
            "covered"
          )]
        ),
        digits = 3,
        row.names = FALSE
      )
      cat(sprintf(
        "  fit_time %.1f s, %d warnings\n",
        res$fit_time[1],
        res$n_warnings[1]
      ))
    }
  }
}

# ---- ic_par closed-form cross-check vs getFitEsts ------------------------------
cat("\n== ic_par closed-form vs getFitEsts cross-check ==\n")
cell <- CELLS[CELLS$cell_id == "core-const-ph-random-r1.5-n300-m10", ]
dat <- generate_data(cell, 42L)
fit <- icenReg::ic_par(
  Surv(L, R, type = "interval2") ~ x_num,
  data = dat$icd,
  model = "ph",
  dist = "weibull"
)
cf <- fit$coefficients
off <- as.numeric(fit$covarOffset)
S_manual <- exp(
  -exp(exp(cf[1]) * (log(TEVAL) - cf[2])) * exp(cf["x_num"] * (0 - off))
)
S_pkg <- 1 -
  icenReg::getFitEsts(fit, newdata = data.frame(x_num = 0), q = TEVAL)
cat("manual:", round(S_manual, 4), "\n")
cat("pkg:   ", round(S_pkg, 4), "\n")
check(
  max(abs(S_manual - S_pkg)) < 1e-6,
  "ic_par closed form matches getFitEsts"
)

cat(sprintf(
  "\n%s\n",
  if (ok_all) "ALL CHECKS PASSED" else "SOME CHECKS FAILED"
))
quit(status = as.integer(!ok_all))
