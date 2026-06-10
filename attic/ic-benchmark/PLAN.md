# Monte Carlo Benchmark: Interval-Censored PAMMs (pamm_ic) — Validation & Comparison

## Context

The `claude/keen-dirac-ZfYvO` branch of pammtools adds interval-censored (IC) survival
modeling via multiple imputation: `pamm_ic()` / `pamm_ic_cr()` (R/pamm-ic.R), conditional-hazard
imputation (R/impute-ic.R), and Rubin-pooled inference (`add_hazard/add_surv_prob/add_cumu_hazard/add_cif.pamm_ic`,
R/pool-ic.R). A 200-rep prototype study exists at `attic/ic-simulation-study.R`.
Before merging/publishing we need a rigorous ADEMP-structured benchmark that
(A1) validates coverage/bias of the MI workflow and (A2) compares it against external IC
methods (icenReg, Turnbull NPMLE), including misspecified-by-design cells.
Production runs on LRZ CoolMUC-4 via SLURM job arrays (R 4.3.3, serial_std partition).

This plan was reviewed by the council of bots (Codex, Gemini, Claude); all consensus
fixes are incorporated below (marked ⊕ where they changed the original design).

**Blocker noted**: LRZ SSH tunnel currently down — user must run `ssh lrz` (password + 2FA)
and keep the session open before Phase 3 (LRZ pilot). All earlier phases are local.

## Locked design decisions (user sign-off)

- Scope: single-event factorial with external comparators + small CR arm (internal comparators only).
- n_sim = **500 reps/cell** (MC SE for 95% coverage ≈ 1.0%). Gemini argued for ≥1900;
  user chose 500. ⊕ Mitigation: headline method comparisons use **paired differences**
  (same data per rep) with paired MC SEs — far more precise than per-method coverage.
- **Binary x ∈ {0,1}** (p=0.5) in the single-event factorial; continuous x only in the
  smooth-f(x) spoke and CR arm.
- `ic_sp` bootstrap CIs in a **4-cell subset** (⊕ corrected from "6": {const,peaked} ×
  PH × rate {1.5, 0.3} × n=300; ≈200k extra ic_sp fits — pilot timing decides
  bs_samples ∈ {50, 100} or shrinks further). `ic_par` = cheap external coverage
  comparator everywhere applicable. Turnbull (`icenReg::ic_np` stratified by x) =
  point-accuracy only, scored at midpoint of NPMLE equivalence class, half-gap +
  fraction-of-TEVAL-in-wide-gaps (>0.2) recorded. ⊕ `survival::survfit` DOES accept
  type="interval2" (verified by council) — use as independent Turnbull cross-check in P0.
- Location: `attic/ic-benchmark/`; add `^attic$` to `.Rbuildignore`. Results gitignored,
  fetched from LRZ via scp/rsync.

## ADEMP (→ attic/ic-benchmark/README.md, locked at Gate R2)

- **Aims**: A1 validation of pamm_ic MI inference; A2 comparison vs ic_par/ic_sp/Turnbull
  incl. misspecified-by-design cells. ⊕ The *headline* misspecification exhibit is
  ic_par's S(t|x) coverage under the peaked baseline (Weibull can't bend) — S(t|x) is
  well-defined for every method. PH-comparator behavior in TV cells is reported as
  estimand mismatch (see Estimands), never as headline "coverage failure".
- **DGPs**: `sim_pexp()` + `add_inspections()`; horizon 10, FINE = seq(0,10,0.05),
  analysis grid CUT = seq(0,10,0.5) frozen everywhere, TEVAL = c(1,2,3,5,7) ⊂ CUT.
- **Estimands** (⊕ all aligned to the piecewise-exponential DGP — see Truth below):
  - S(t|x) at TEVAL × x∈{0,1} (exactly aligned: PAMM survival cumulates h_j·intlen).
  - Hazard: **interval-average** h̄_j = (H(t_j)−H(t_{j−1}))/0.5 on CUT intervals
    containing TEVAL (⊕ was: smooth h(t) at endpoints — council showed endpoint vs
    interval-average differs by up to ~25% on the peak flanks and would corrupt
    hazard coverage for ALL methods). ic_par scored the same way via its closed-form
    cumulative hazard difference.
  - logHR(t): interval-average log h̄(t|1) − log h̄(t|0); equals β under PH. ⊕ In TV
    cells, PH-only methods (ic_par/ic_sp) are scored as **projection error vs β(t)**
    plus bias against the **least-false constant logHR** (computed once per TV cell
    from one large-n fit at config time) — not as CI coverage.
  - β scalar (PH cells only); f(x) on x-grid, ⊕ estimate AND truth both centered to
    mean zero over the evaluation grid before scoring (spoke S); CIF_k(TEVAL) at
    ⊕ x = 0 (CR arm, as in prototype).
- **Methods & model formulas** (⊕ per-cell formulas now explicit — pamm_ic's default
  `~ s(tend) + x` is a PH model and cannot recover β(t)):
  - PH cells: `ped_status ~ s(tend, k = 10) + x`.
  - TV cells: `ped_status ~ s(tend, k = 10) + x + s(tend, by = x_num, k = 10)`
    for pamm_ic, midpoint AND oracle (x_num = numeric 0/1; by-term carries its own k).
  - Spoke S: `ped_status ~ s(tend, k = 10) + s(x, k = 10)`.
  - CR arm: prototype's `ped_status ~ s(tend, by = cause) + cause + cause:x`.
  - External: ic_par (Weibull PH), ic_sp (PH), ic_np — single-event linear-effect cells only.
- **Performance measures**: bias, empirical SE, RMSE, coverage, bias-eliminated coverage,
  mean CI width, ⊕ avSE/empSE ratio (model-based SE column added to schema — the most
  diagnostic quantity for explaining coverage failures), failure/non-convergence rate
  (incl. ⊕ icenReg convergence flags, mgcv warnings), fit_time. ⊕ MC SEs for term-level
  (grid-averaged) coverage computed by **clustering at the rep level**
  (sd of rep-level means / √500), never pooled-Bernoulli. Paired method-difference
  MC SEs for headline contrasts.

### Truth computation (⊕ revised: exact, not integrate())

`sim_pexp()` draws from a **left-endpoint piecewise-constant hazard on FINE** — not from
the smooth formula. Therefore: config.R defines each DGP as a *rate vector on FINE*
(evaluated from the log-hazard function at left endpoints), and BOTH the `sim_pexp`
formula and all truths derive from it. S(t), interval-average h̄_j, logHR(t), least-false
projections, and CIF_k (CR) are computed **exactly** via cumulative sums of the FINE
rates — no `integrate()`, no Riemann drift, truth and simulator literally cannot diverge.
P0 still cross-checks against a 10⁶-draw KM.

### Factorial: core + spokes (⊕ 39 cells × 500 reps = 19,500 reps)

| Arm | Factors | Cells |
|---|---|---|
| Core (mechanism=random) | baseline {const: −1.5; peaked: −3.2+1.6·exp(−((t−3)/1)²/2)} × effect {ph: +0.6·x; tv: +1.2·exp(−0.3t)·x} × rate {1.5, 0.6, 0.3} × n {300, 1000} | 24 |
| Spoke M (mechanism) | {fixed, mixed} × {dense, sparse schedule}; peaked, PH, n=300 | 4 |
| Spoke R (rising baseline) | −3.7+0.12t, PH, n=300, 3 rates | 3 |
| Spoke S (smooth f(x)) | peaked + 0.8·sin(πx), x~U(−1,1); internal methods only; n=300, rate {1.5, 0.3} | 2 |
| ⊕ Spoke I (imputation count) | m ∈ {5, 20} at peaked, PH, n=300, rate 0.6 — informs the package default m | 2 |
| CR arm | prototype truths (cause1 peaked B1=0.7, cause2 rising B2=−0.5), x~U(−1,1), n=450, 3 rates; MI/midpoint/oracle | 3 |
| ⊕ CR-unknown | one CR cell (rate 0.6) with 30% of event causes masked — exercises `pamm_ic_cr`'s cause-imputation branch (R/pamm-ic.R:240-247), otherwise completely unvalidated | 1 |
| (subset flag) | ic_sp bootstrap in 4 core cells (above) | — |

Known acknowledged gaps (README): mechanism × n and mechanism × effect interactions
untested; spoke conclusions are conditional on the (peaked, PH, n=300) reference.

Compute estimate: ~250 CPU-h (plan 400 with margin); replaced by measured pilot timings.

## File layout (attic/ic-benchmark/)

```
README.md             ADEMP + run instructions + design-decision log
config.R              factor levels, cells table (incl. per-cell model_formula,
                      FINE rate vectors, exact truths, least-false projections),
                      task_table (cell × rep) with seeds; pins RNGkind; vendors
                      rubin_scalar from the prototype
dgp.R                 rate-vector builders, exact truth functions, generate_data(cell, seed)
methods.R             fit_<method>() + extract_<method>() returning standardized
                      long tibble; get_applicable_methods(cell); bespoke pooled
                      logHR(t) helper (see below)
metrics.R             score_rep(); rep-clustered + paired MC-SE helpers
run-task.R            CLI: Rscript run-task.R <task_id> <n_cores>
run-local.R           smoke/pilot driver (local mclapply)
phase0-known-answer.R known-answer validation
debug-single-fit.R    one (cell, rep, method) verbose, replayable from stored seeds
aggregate.R           bind raw, defensive dedup (assert identical est/lower/upper/
                      truth on duplicate keys, exclude fit_time), completeness check
                      vs task_table, pointwise + term-level summaries → summary.{rds,csv}
report.qmd            figures + computed-not-speculative prose (reads summary.rds only)
slurm/install-lrz.sh  module load r/4.3.3-gcc13-mkl; git pull branch; R CMD INSTALL;
                      install icenReg + deps; ⊕ verify serial_std grants
                      --cpus-per-task=16 on one node (scontrol show part)
slurm/production.sbatch  --clusters=serial --partition=serial_std --array=1-39
                      --cpus-per-task=16 --mem=24G --export=NONE --get-user-env;
                      OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE;
                      ⊕ per-task --time from pilot s/rep × 500/16 × 2 (load
                      imbalance spans ~an order of magnitude across cells)
slurm/pilot.sbatch    IC_PILOT_REPS=20, flat --time=00:40:00
results/, figures/    gitignored
```

### Execution architecture (⊕ revised after council)

- One array task = one cell (40 tasks). Inside: ⊕ **wave-batched** `mclapply`
  (chunks of ~32 reps, `mc.preschedule = FALSE`); the **parent** writes
  `results/raw/task_<id>.rds` after each wave via atomic tempfile+rename — forked
  children never write shared files. Resume keys on completed rep ids in the file.
  `try-error` elements from mclapply explicitly checked and converted to failure rows.
- ⊕ **Seed substreams**: task_table stores `seed_data` per rep; each stochastic method
  gets a deterministic sub-seed (`seed_data + method_index · 1e7`-style, set before the
  method call). Adding/removing a method never changes other methods' results;
  completed cells stay valid under design edits.
- Per-method tryCatch (errors → row with error_msg, never dropped); warnings captured.
  Replicate functions return ONLY metric rows — model objects discarded before return
  (memory guard for 16 workers × m=10 fits).
- Row schema: cell_id, rep, seed_data, method, estimand, t, x, truth, est, **se**,
  lower, upper, covered, width, err, misspecified, error_msg, n_warnings, conv_flag,
  fit_time, mean_int, prop_right_cens, ⊕ prop_left_cens, n_inspections_median.
  ~20k reps × ~75 rows ≈ 1.5M rows, tens of MB.

## Method extraction specifics

- ⊕ All simulation-based CIs use **nsim = 500 draws** (pamm_ic pooled draws, midpoint/
  oracle `ci_type="sim"`, ic_par pushed-through normal draws) — equalized so coverage
  differences aren't tail-order-statistic artifacts; with m=10 that's 50 draws/fit
  before Rubin inflation. Pilot includes one nsim ∈ {200, 500, 1000} sensitivity cell.
- Survival/CIF prediction grids ⊕ **must be `group_by(x)`** (resp. `group_by(cause, x)`)
  before `add_surv_prob`/`add_cif` — ungrouped grids silently accumulate hazard across
  strata (R/add-functions.R:201). P1 asserts S(t|0/1) monotone within stratum.
- Bespoke pooled logHR(t) helper: per imputation fit, paired draws
  `(X1 − X0) %*% beta_draw` from the SAME coefficient draw (never differences of
  independently simulated quantities), then `rubin_inflate_qoi_draws`-style inflation
  on the logHR scale, mirroring R/pool-ic.R:58-71. ⊕ Explicit Gate-R1 checklist item.
- β CIs: pamm_ic via rubin_scalar with ⊕ **Barnard–Rubin small-sample df** (m=10 with
  high FMI makes classic-Rubin-t vs z a material construction difference); df
  distribution reported. midpoint/oracle/ic_par: z-intervals; construction difference
  noted wherever β coverage is compared.
- ic_par: closed-form Weibull h̄_j (via cumulative hazard differences) and S;
  cross-check vs getFitEsts/survCIs in P1.
- ic_sp: getFitEsts S points; β point always; bootstrap CIs (bsMat) in the 4-cell subset.
- Turnbull ic_np: getSCurves midpoint convention (above).
- ⊕ Clipping diagnostic: fraction of pooled draws clipped at [0,1]/0 boundaries
  (R/pool-ic.R:124-132) recorded — matters for S(7|x) ≈ 0.06 tails.

## Phases & review gates

| Phase | Where | Content / exit criterion |
|---|---|---|
| P0 known-answer | local | oracle PAMM, const+PH, n=1000, 200 reps → β & S(t) coverage within ⊕ ±3·MC SE of 0.95 (≈[90.4, 99.6]%, diagnostic not brittle gate); pamm_ic with near-exact intervals (rate=20) ≈ oracle; exact truth vs 10⁶-draw KM (max diff < 0.005); ⊕ survfit-vs-ic_np Turnbull cross-check; icenReg API checks |
| P1 single-fit debug | local | every method × extractor on one rep; schema conformance; ic_par CI cross-check; ⊕ grouped-grid monotonicity assertion; logHR paired-draw sanity (PH cell: logHR(t) flat ≈ β̂) |
| **Gate R1** | — | council-of-bots review of config.R/dgp.R/methods.R/metrics.R/run-task.R; checklist: paired logHR draws, Rubin inflation scale, exact-truth derivation, seed substreams, wave-write architecture; fix HIGH/CRITICAL; re-run P0 |
| P2 smoke | local | 2 reps × 5 cheap cells (one per arm incl. CR-unknown), m=3, nsim=50 via IC_SMOKE=1; full path incl. aggregate.R + report render; exercise resume (kill + rerun) |
| P3 pilot + timing | **LRZ** | pilot.sbatch, 20 reps × 39 cells; measured s/rep → per-task --time + final CPU-h; sanity: oracle ≈ nominal everywhere, midpoint under-covers sparse cells, failure rate <5%/cell; finalize ic_sp bs_samples/subset; nsim sensitivity check |
| **Gate R2** | — | council-of-bots on pilot results; ADEMP locked; ⊕ pre-commit failure-rate rule (flag cells >2% failures — survivor-bias guard); post-pilot changes → decision log |
| P4 production | **LRZ** | production.sbatch array 1-39, 500 reps; resubmit failed indices only; fetch results/raw |
| P5 aggregate + report | local | aggregate.R + report.qmd (all prose computed) |
| **Gate R3** | — | council-of-bots on summary + figures before reporting |

LRZ prerequisites (before P3): user re-establishes SSH tunnel; install-lrz.sh once;
git push/pull sync for code; results back via scp.

## Report figures

1. (headline) Coverage of S(t|x) & h̄(t|x) vs mean inspection-interval length, faceted
   estimand × baseline/effect, colored by method, rep-clustered MC-SE pointranges, 0.95 line.
2. β / logHR(t) coverage incl. ic_par & ic_sp-boot; ⊕ paired MI−midpoint coverage
   differences with paired MC SEs as inset/panel.
3. RMSE + RRMSE-vs-oracle of S(t) points incl. Turnbull (only Turnbull figure;
   wide-gap fraction annotated).
4. CI width + ⊕ avSE/empSE ratio (mechanism: midpoint too narrow, MI inflated correctly).
5. TV cells: estimated vs true β(t); ⊕ ic_par least-false constant overlay
   (decomposes its error into estimand mismatch vs estimation noise).
6. CR arm coverage/RMSE incl. CR-unknown cell; spoke-S f(x) curves; spoke-I m-sensitivity.
Optional: rsimsum zip plots; nested-loop plot over core factorial.
Tables: term-level summaries (MC SEs), failure rates, median fit_times, FMI summaries.

## Known risks (tracked in README)

- ic_sp bootstrap cost (~200k extra NPMLE fits in subset) — pilot decides bs_samples
  (100→50) or drops ic_sp coverage (ic_par remains).
- CR add_cif cost — fallback nsim 500→200 for CIF CIs only.
- icenReg source compile on LRZ (C++/gcc13); icenReg may spawn own threads —
  verify single-threaded behavior during install (council flag).
- pammtools deps (mvtnorm, scam, pec, Formula, ...) in LRZ 4.3 user lib.
- Per-cell walltime heterogeneity — set from pilot, not flat.

## Verification

- P0 known-answer thresholds are the infrastructure test (±3·MC SE diagnostic gates).
- P2 smoke must produce a rendered report from real (tiny) output end-to-end, incl. resume.
- P3 pilot must reproduce prototype's qualitative findings (midpoint under-coverage
  worsening with sparsity; MI tracking oracle) before production launches.
- Gates R1–R3 (council-of-bots) at code, pilot, and results stages.
