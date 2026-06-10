# ===========================================================================
# IC benchmark: configuration (single source of truth for the design)
# ---------------------------------------------------------------------------
# Defines all constants, the DGP cell table, the method registry, and the
# task table (cell x rep) with upfront seeds. Sourced by every other script.
# See README.md for the ADEMP design document.
# ===========================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

# RNG pinned so stored seeds reproduce across machines / R builds (R >= 4.1)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ---- constants ------------------------------------------------------------
BASE_SEED <- 20260610L
HORIZON <- 10
FINE_STEP <- 0.05
FINE <- seq(0, HORIZON, by = FINE_STEP) # simulation grid (sim_pexp cut)
CUT <- seq(0, HORIZON, by = 0.5) # analysis grid, frozen everywhere
TEVAL <- c(1, 2, 3, 5, 7) # evaluation times (subset of CUT)
stopifnot(all(TEVAL %in% CUT), all(CUT %in% FINE))

M_DEFAULT <- 10L # imputations (spoke I overrides per cell)
NSIM <- 500L # posterior draws for ALL simulation-based CIs (equalized)
NREP_MAX <- 500L # production reps per cell; seeds always generated for this
ALPHA <- 0.05
K_TEND <- 10L # basis dimension for s(tend) terms
X_GRID <- seq(-0.9, 0.9, by = 0.2) # evaluation grid for smooth f(x), spoke S
# ic_sp bootstrap draws (subset cells only): 500 so percentile-CI quantile
# noise is commensurable with the 500 posterior draws used everywhere else;
# pilot timing may lower this (recorded in the decision log if so)
BS_SAMPLES_ICSP <- 500L

# runtime overrides (smoke / pilot)
IC_SMOKE <- nzchar(Sys.getenv("IC_SMOKE"))
N_REP <- if (IC_SMOKE) {
  2L
} else if (nzchar(Sys.getenv("IC_PILOT_REPS"))) {
  as.integer(Sys.getenv("IC_PILOT_REPS"))
} else {
  NREP_MAX
}
M_RUN <- if (IC_SMOKE) 3L else NA_integer_ # NA = use cell$m
NSIM_RUN <- if (IC_SMOKE) 50L else NSIM

# ---- truth parameters (shared with dgp.R) ----------------------------------
BETA_PH <- 0.6
B1_CR <- 0.7
B2_CR <- -0.5

# ---- cells -----------------------------------------------------------------
# model_type: "ph" | "tv" | "smooth" | "cr"
# Inspection schedules for spoke M (mean gap 2/3 = "dense", 10/3 = "sparse")
SCHED_DENSE <- seq(2 / 3, HORIZON, by = 2 / 3)
SCHED_SPARSE <- seq(10 / 3, HORIZON, by = 10 / 3)

cells_core <- expand_grid(
  baseline = c("const", "peaked"),
  effect = c("ph", "tv"),
  rate = c(1.5, 0.6, 0.3),
  n = c(300L, 1000L)
) |>
  mutate(
    arm = "core",
    mechanism = "random",
    schedule = list(NULL),
    model_type = effect,
    m = M_DEFAULT,
    icsp_boot = baseline %in%
      c("const", "peaked") &
      effect == "ph" &
      rate %in% c(1.5, 0.3) &
      n == 300L,
    cause_mask = 0
  )

cells_mech <- expand_grid(
  mechanism = c("fixed", "mixed"),
  sched_lab = c("dense", "sparse")
) |>
  mutate(
    baseline = "peaked",
    effect = "ph",
    n = 300L,
    # for "mixed", add_inspections jitters by +U(0, 1/rate): set rate = 1/gap
    # AND start the schedule at 0 so the offset phase-randomizes the grid
    # within the first gap (a schedule starting at `gap` would push the first
    # inspection into (gap, 2*gap] -- R2 finding)
    rate = ifelse(sched_lab == "dense", 1.5, 0.3),
    schedule = map2(
      sched_lab,
      mechanism,
      \(s, mech) {
        gap <- if (s == "dense") 2 / 3 else 10 / 3
        if (mech == "mixed") {
          seq(0, HORIZON, by = gap)
        } else {
          seq(gap, HORIZON, by = gap)
        }
      }
    ),
    arm = "mech",
    model_type = "ph",
    m = M_DEFAULT,
    icsp_boot = FALSE,
    cause_mask = 0
  ) |>
  select(-sched_lab)

cells_rising <- tibble(rate = c(1.5, 0.6, 0.3)) |>
  mutate(
    baseline = "rising",
    effect = "ph",
    n = 300L,
    arm = "rising",
    mechanism = "random",
    schedule = list(NULL),
    model_type = "ph",
    m = M_DEFAULT,
    icsp_boot = FALSE,
    cause_mask = 0
  )

cells_smooth <- tibble(rate = c(1.5, 0.3)) |>
  mutate(
    baseline = "peaked",
    effect = "smooth",
    n = 300L,
    arm = "smooth",
    mechanism = "random",
    schedule = list(NULL),
    model_type = "smooth",
    m = M_DEFAULT,
    icsp_boot = FALSE,
    cause_mask = 0
  )

cells_m <- tibble(m = c(5L, 20L)) |>
  mutate(
    baseline = "peaked",
    effect = "ph",
    rate = 0.6,
    n = 300L,
    arm = "mimp",
    mechanism = "random",
    schedule = list(NULL),
    model_type = "ph",
    icsp_boot = FALSE,
    cause_mask = 0
  )

cells_cr <- tibble(rate = c(1.5, 0.6, 0.3), cause_mask = 0) |>
  bind_rows(tibble(rate = 0.6, cause_mask = 0.3)) |>
  mutate(
    baseline = "cr",
    effect = "cr",
    n = 450L,
    arm = ifelse(cause_mask > 0, "cr_unknown", "cr"),
    mechanism = "random",
    schedule = list(NULL),
    model_type = "cr",
    m = M_DEFAULT,
    icsp_boot = FALSE
  )

CELLS <- bind_rows(
  cells_core,
  cells_mech,
  cells_rising,
  cells_smooth,
  cells_m,
  cells_cr
) |>
  mutate(
    cell_id = sprintf(
      "%s-%s-%s-%s-r%s-n%d-m%d%s",
      arm,
      baseline,
      effect,
      mechanism,
      format(rate),
      n,
      m,
      ifelse(cause_mask > 0, "-mask", "")
    ),
    task_id = row_number()
  ) |>
  relocate(task_id, cell_id)
stopifnot(nrow(CELLS) == 39L, !anyDuplicated(CELLS$cell_id))

# ---- model formulas (explicit per model_type) ------------------------------
# x is a factor with numeric companion x_num (needed for the TV by-term)
model_formula_for <- function(model_type) {
  switch(
    model_type,
    ph = ped_status ~ s(tend, k = 10) + x,
    tv = ped_status ~ s(tend, k = 10) + x + s(tend, by = x_num, k = 10),
    smooth = ped_status ~ s(tend, k = 10) + s(x, k = 10),
    cr = ped_status ~ s(tend, by = cause, k = 10) + cause + cause:x,
    stop("unknown model_type: ", model_type)
  )
}

# ---- method registry --------------------------------------------------------
# Fixed indices: sub-seeds are derived per (rep, method index), so adding new
# methods (new indices) never changes existing methods' draws.
METHOD_INDEX <- c(
  mi = 1L,
  midpoint = 2L,
  oracle = 3L,
  ic_par = 4L,
  ic_sp = 5L,
  turnbull = 6L
)

get_applicable_methods <- function(cell) {
  switch(
    cell$model_type,
    ph = c("mi", "midpoint", "oracle", "ic_par", "ic_sp", "turnbull"),
    tv = c("mi", "midpoint", "oracle", "ic_par", "ic_sp", "turnbull"),
    smooth = c("mi", "midpoint", "oracle"),
    cr = c("mi", "midpoint", "oracle")
  )
}

# ---- task table + seeds -----------------------------------------------------
# Seeds are ALWAYS generated for NREP_MAX reps so smoke/pilot subsets use the
# same per-rep seeds as production (resumable, extendable).
set.seed(BASE_SEED)
TASK_TABLE <- expand_grid(task_id = CELLS$task_id, rep = seq_len(NREP_MAX)) |>
  mutate(seed_data = sample.int(.Machine$integer.max, n())) |>
  left_join(CELLS |> select(task_id, cell_id), by = "task_id")

# Sub-seeds for stochastic methods: deterministic per (rep, method), derived
# from seed_data without touching the data-generation stream order.
method_subseeds <- function(seed_data, n_methods = 32L) {
  set.seed(seed_data)
  sample.int(.Machine$integer.max, n_methods)
}

# ---- paths ------------------------------------------------------------------
BENCH_DIR <- {
  # robust whether sourced from package root or from attic/ic-benchmark
  cand <- c("attic/ic-benchmark", ".")
  cand[file.exists(file.path(cand, "config.R"))][1]
}
stopifnot(!is.na(BENCH_DIR))
RAW_DIR <- file.path(BENCH_DIR, "results", "raw")
RES_DIR <- file.path(BENCH_DIR, "results")
FIG_DIR <- file.path(BENCH_DIR, "figures")
for (d in c(RAW_DIR, FIG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
