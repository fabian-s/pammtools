# IC Benchmark: Monte Carlo validation & comparison of interval-censored PAMMs

Simulation study for the `pamm_ic()` / `pamm_ic_cr()` multiple-imputation (MI)
workflow. Designed per the ADEMP framework (Morris, White & Crowther 2019).
See `PLAN.md` for the full reviewed implementation plan and
`attic/ic-simulation-study.R` for the earlier 200-rep prototype.

## ADEMP

**Aims.**
A1: validate that `pamm_ic` MI inference attains nominal 95% coverage and low
bias for hazard, survival, covariate effects and time-varying effects under
interval censoring, across inspection sparsity, inspection mechanisms and
baseline shapes.
A2: compare against external IC methods — `icenReg::ic_par` (Weibull PH),
`icenReg::ic_sp` (semiparametric PH), Turnbull NPMLE (`icenReg::ic_np`) —
including **misspecified-by-design** cells. The headline misspecification
exhibit is ic_par's S(t|x) coverage under non-Weibull (peaked) baselines;
PH-only comparators in time-varying (TV) cells are reported as estimand
mismatch, not headline "coverage failure".

**Data-generating mechanisms.** `sim_pexp()` + `add_inspections()`. Horizon 10;
simulation grid FINE = seq(0, 10, 0.05); analysis grid CUT = seq(0, 10, 0.5)
frozen across all cells and imputations; evaluation times TEVAL = 1, 2, 3, 5, 7.
39 cells x 500 reps (see Design below).

**Estimands** (all aligned exactly to the piecewise-exponential DGP):

- S(t|x) at TEVAL x {0,1}
- hazard as the **interval-average** hbar_j = (H(t_j) - H(t_{j-1})) / 0.5 on the
  CUT interval ending at each TEVAL (endpoint-vs-interval-average mismatch would
  corrupt hazard coverage for all methods; see PLAN.md)
- logHR(t) = log hbar_j(1) - log hbar_j(0); equals beta under PH
- beta (PH cells); centered smooth f(x) on a grid (spoke S); CIF_k(t) at x=0 (CR arm)

**Methods.** `pamm_ic` (MI, proper, m per cell), midpoint-imputation PAMM,
oracle PAMM (exact times), ic_par, ic_sp, Turnbull. Applicability: external
methods only in single-event linear-effect cells; CR and smooth arms are
internal-only. Turnbull is point-accuracy only (NPMLE has no standard CIs),
scored at the midpoint of its equivalence class with the half-gap recorded as
`ambiguity`. ic_sp survival is point-accuracy only; its beta gets bootstrap
percentile CIs (bs_samples = 100) in 4 designated cells.

**Performance measures.** Bias, empirical SE, RMSE, coverage, bias-eliminated
coverage, mean CI width, avSE/empSE ratio, failure/non-convergence rate,
fit time — with Monte Carlo SEs. Term-level (grid-averaged) coverage MC SEs are
**clustered at the rep level**; headline MI-vs-midpoint contrasts use **paired
differences** (same data per rep). n_sim = 500/cell: MC SE for 95% coverage
~ 0.97%; paired differences are substantially more precise.

## Design: core factorial + spokes (39 cells)

| Arm | Factors | Cells |
|---|---|---|
| core | baseline {const −1.5; peaked −3.2+1.6·exp(−((t−3)/1)²/2)} × effect {ph +0.6x; tv +1.2·exp(−0.3t)·x} × rate {1.5, 0.6, 0.3} × n {300, 1000}; mechanism random | 24 |
| mech | mechanism {fixed, mixed} × schedule {dense (gap 2/3), sparse (gap 10/3)}; peaked, ph, n=300 | 4 |
| rising | baseline −3.7+0.12t; ph, n=300, 3 rates | 3 |
| smooth | peaked + 0.8·sin(πx), x~U(−1,1); internal methods; n=300, rate {1.5, 0.3} | 2 |
| mimp | m ∈ {5, 20}; peaked, ph, n=300, rate 0.6 | 2 |
| cr | cause1 peaked (β1=0.7), cause2 rising (β2=−0.5); x~U(−1,1), n=450, 3 rates | 3 |
| cr_unknown | as cr (rate 0.6) with 30% of event causes masked (exercises cause imputation) | 1 |

x is binary (p = 0.5) in single-event cells so stratified Turnbull is directly
comparable. Acknowledged gaps: mechanism × n and mechanism × effect
interactions untested; spoke conclusions are conditional on the
(peaked, ph, n=300) reference cell.

## Truth computation

`sim_pexp()` draws from a left-endpoint piecewise-constant hazard on FINE — not
from the smooth formula. All truths (S, interval-average hazard, logHR, CIF)
are computed **exactly** from those same FINE rate vectors via cumulative sums
(`dgp.R`), so truth and simulator derive from one definition. Phase 0
cross-checks against large-sample Kaplan–Meier.

## How to run

```bash
# Phase 0 known-answer validation (local, ~15-30 min)
Rscript attic/ic-benchmark/phase0-known-answer.R

# P1 single-fit debugging (local, minutes)
Rscript attic/ic-benchmark/debug-single-fit.R

# P2 smoke test (local, minutes): 2 reps x 5 cheap cells, m=3, nsim=50
IC_SMOKE=1 Rscript attic/ic-benchmark/run-local.R
Rscript attic/ic-benchmark/aggregate.R

# P3 pilot on LRZ (after bash slurm/install-lrz.sh on the login node)
sbatch attic/ic-benchmark/slurm/pilot.sbatch

# P4 production on LRZ
sbatch attic/ic-benchmark/slurm/production.sbatch
# resubmit individual failed cells: sbatch --array=<ids> [--time=...] ...

# P5 aggregation + report (local, after fetching results/raw from LRZ)
Rscript attic/ic-benchmark/aggregate.R
quarto render attic/ic-benchmark/report.qmd
```

Seeds: one base seed generates per-rep `seed_data` for all 500 reps upfront
(`config.R`); each method gets a deterministic sub-seed derived from `seed_data`,
so adding/removing methods never changes other methods' results, and smoke/
pilot subsets use the same seeds as production. RNGkind is pinned.

## Design-decision log

| Date | Decision |
|---|---|
| 2026-06-10 | Initial design approved (PLAN.md), incl. council-review fixes: interval-average hazard estimand; exact PC-grid truths; nsim equalized at 500; per-(rep,method) sub-seeds; rep-clustered MC SEs; ic_sp bootstrap restricted to 4 cells; Turnbull point-only with midpoint-of-gap convention; PH-comparator results in TV cells reported as estimand mismatch. |
| 2026-06-10 | ic_par closed forms must un-center covariates via `fit$covarOffset` (icenReg centers internally; verified against `getFitEsts`). |
| 2026-06-10 | ic_sp S(t|x): point accuracy only (bsMat covers coefficients, not the NPMLE baseline, so honest bootstrap S-CIs are not available). ic_par is the external S-coverage comparator. |
| 2026-06-10 | Least-false constant logHR for TV cells: computed from one large-n exact-time fit (coxph / survreg) per cell — exact-time projection convention (cheap, censoring-free); used as reporting overlay only. |
| 2026-06-10 | Gate R1 fixes: Turnbull upper/lower envelope assignment corrected (ambiguity was sign-flipped); ic_sp bootstrap raised to bs_samples=500 for draw-count commensurability with the 500 posterior draws elsewhere; transient (rep/fork-level) failures retried on resume; icenReg loaded in the runner parent. Note: MI beta CIs use Barnard–Rubin t; midpoint/oracle/ic_par use Wald z — a deliberate construction difference, stated wherever beta coverage is compared. |

ADEMP is locked at Gate R2 (pre-production); any later change goes in this log.
