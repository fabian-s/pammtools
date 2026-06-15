#!/usr/bin/env Rscript
# ===========================================================================
# IC benchmark: least-false constant logHR projections for the TV cells
# ---------------------------------------------------------------------------
# One large-n exact-time fit per distinct TV DGP (projection depends on the
# baseline + effect + admin censoring, not on inspection rate or n), written
# to results/least-false.rds for the report overlay (Fig 5).
# Usage: Rscript attic/ic-benchmark/compute-least-false.R
# ===========================================================================

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

tv_cells <- CELLS |> filter(effect == "tv")
distinct_dgps <- tv_cells |> distinct(baseline)

lf <- purrr::map_dfr(distinct_dgps$baseline, function(bl) {
  cell <- tv_cells |> filter(baseline == bl) |> slice(1)
  out <- compute_least_false(cell, n = 1e5, seed = 4242L)
  out$baseline <- bl
  out
})
# expand to all TV cell_ids
lf_all <- tv_cells |>
  select(cell_id, baseline) |>
  left_join(lf |> select(-cell_id), by = "baseline")

saveRDS(lf_all, file.path(RES_DIR, "least-false.rds"))
print(as.data.frame(lf_all), digits = 3)
cat("wrote", file.path(RES_DIR, "least-false.rds"), "\n")
