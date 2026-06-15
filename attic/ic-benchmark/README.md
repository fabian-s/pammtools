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

**Fixed-schedule identifiability (R2 annotation).** In the fixed-schedule
cells every subject shares the same inspection grid, so S(t) is
nonparametrically identified ONLY at the schedule points; at other TEVAL
points all methods score model-based interpolation, and hazard estimands can
be degenerate (vacuous CIs, extreme estimates) on the sparse 3-point grid.
These cells are robustness exhibits and are never headlined alongside the
random-inspection cells.

**Iterated-MI extension (scope note).** `pamm_ic` draws all m imputations
from the single midpoint-initialized fit (one-step MI). The pilot suggests
this leaves early-time bias under sparse inspection. Whether one
refit-and-reimpute iteration removes that bias is OUT OF SCOPE for this
study; thanks to the fixed per-(rep, method-index) seed substreams it can be
added later as a clean extension over the same seeds without invalidating
production results.

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
| 2026-06-10 | Gate R2 (pilot, 20 reps × 39 cells, zero failures): (1) term-level avSE/empSE recomputed within (t,x) grid points and median-summarized — pooling errors across the grid mixed bias heterogeneity into empSE (the pilot's apparent se_ratio 0.44 was an aggregation artifact; pointwise ratios ≈ 1). (2) Mixed-mechanism schedules now start at 0 so the +U(0, gap) offset phase-randomizes within the first gap; stale pilot results for the two mixed cells (tasks 27, 28) deleted before production. (3) PLAN's `conv_flag`/`n_inspections_median` columns dropped from the schema (mgcv warnings + failure rates cover non-convergence; inspection density is summarized by mean_int and censoring proportions); PLAN's clipping-fraction diagnostic realized as the boundary_lo/hi proxy. (4) Flat 4 h walltime for all 39 tasks (slowest projected ≈ 45 min) instead of per-task limits. (5) Sparse-cell MI under-coverage diagnosed as honest one-step-MI point-estimate bias (bias-eliminated coverage nominal; oracle nominal; Turnbull recovers S(1) on the same data) — a study finding, not an implementation issue. (6) Spoke-I (m) and MI≈oracle equivalence statements are precision-limited at 500 reps (±1–2% coverage) — reported as qualitative. |

| 2026-06-10 | P0 outcome: ALL statistical gates passed (oracle & near-exact-MI coverage 0.937–0.942 for surv/hazard/beta, exact-truth vs KM ≤ 0.0016, Turnbull cross-check, icenReg API invariants). The 98%-completion criterion failed locally (191/200) due to fork OOM kills on the 15 GB workstation — no statistical failures; rerun submitted on LRZ (job 5256433) for the record. Production proceeds on the statistical verdict. |

| 2026-06-10 | nsim sensitivity (LRZ job 5256447, 20 reps, peaked/ph/r0.6/n300): coverage identical across nsim ∈ {200, 500, 1000} (max diff 0.015, ≪ MC SE) for all methods × estimands — production nsim = 500 confirmed. **ADEMP FROZEN; production launched** (500 reps × 39 cells; reps 1–20 resume from the pilot in 37 cells, mixed-mechanism cells re-run fresh under the corrected schedule). |

| 2026-06-11 | **Gate R3 — DGP defect found; corrective amendment + production re-run.** The R3 review (Claude leg, independently confirmed by Gemini and a second adversarial Claude leg; Codex timed out twice) found an informative-censoring bug in `add_inspections` (R/sim-ic.R): survivors were recorded with their exact admin-censoring time `(10, Inf)` while events after the last inspection got `(L_last, Inf)` — a hybrid coarsening convention violating CAR. Records `(L, Inf)` with `L < 10` actually imply the event lies in `(L, 10]`, but all IC likelihoods score them as `S(L)`; hazard underestimated, S(t) overestimated at late t. Smoking gun: correctly-specified ic_par collapsing to 0.37 S(7) coverage in const/r0.3/n1000. Affected: ~1.5/4.6/12.3% of records at rates 1.5/0.6/0.3; 37/39 cells contaminated (all random + both mixed cells; the two fixed-schedule cells are immune since their schedules end at exactly 10). The Gate-R2 diagnosis "sparse-cell MI under-coverage = honest one-step-MI bias" is amended: TRUE at early t (S(1) bias is genuine one-step-MI bias, inherited from midpoint initialization), but late-t under-coverage and the sparse-cell beta-bias gradients were dominated by this artifact, shared by all non-oracle methods. **Fix:** `add_inspections` gains `terminal_exam = TRUE` (default): every grid additionally contains `max_time`, so undetected events get `R = max_time`; `terminal_exam = FALSE` implements the other coherent convention (survivors censored at their last inspection). The fix is deterministic (no extra RNG draws) — verified bit-identical exact times, grids and L bounds vs the old code; only `R` changes on exactly the undetected events. New: regression test (no event row with `Inf` upper bound below `max_time`), NA-status assertion, `prop_terminal_only` meta diagnostic, and P0 gate (f): well-specified ic_par at const/r0.3 must be nominal at all t (the gate (b)/(c) missed, being near-exact-interval-based). **No P0 re-run pre-production** (user decision): the fix was triple-verified empirically by the council legs and the stream-identity check covers the pipeline; gate (f) guards future changes. Pre-fix results archived as `results-prebugfix/`; full 39-cell production re-run on LRZ. |

| 2026-06-11 | **Gate R3 FINAL: corrected production signed off** (Codex + Gemini + Claude legs, unanimous; no further re-runs). Verified on the corrected run (39 cells × 500 reps, 0 failures): artifact gone exactly where predicted (ic_par const/r0.3/n1000 S(7) coverage 0.372→0.948; beta gradients collapse; oracle bit-identical pre/post, confirming stream identity); `prop_terminal_only` matches predicted contamination (1.3/4.2/12.5% per rate). Post-fix MI/midpoint coverage DROPS in sparse const cells are bias movement, NOT narrower CIs (widths −2–5%; be_coverage nominal throughout): the pre-fix midpoint mid/late-t accuracy was an accidental cancellation with the artifact; on the corrected DGP midpoint's (L+10)/2 imputation of terminal-bracketed events is measured honestly (mechanism verified against a 6e5-subject imputed-data KM limit, matches all 10 (t,x) points). Report fixes applied per council: DGP section rewritten to the terminal-exam convention; computed sparse-cell headline prose (BE coverage + n-gradient + Turnbull S(1) contrast); peaked-baseline hazard-estimand caveat with oracle smoothing ceiling (0.907–0.925); warning-rate disclosure (fixed-sparse MI cell: 93.6% step-failure warnings); MI-vs-midpoint dominance scoped to S/hazard coverage (beta/logHR: 17/36 better, 0 worse; logHR RMSE in TV n300 cells up to 17% worse); nsim-sensitivity sentence anchored to the decision log. |

| 2026-06-11 | **Iterated-MI extension study** (the extension scoped in the one-step note above, run over the same seeds). `pamm_ic` gains `iter` (chained refit-and-reimpute per imputation chain; `iter = 1` = production one-step MI, draw-for-draw). Pilot (10 cells × 50 reps, paired): S(1\|0) bias −35%, term-level surv coverage up in 10/10 cells, no harm anywhere, ~1.6× cost. iter-sensitivity (4 worst cells, iter ∈ {1,2,3,5}): bias decays ~geometrically (worst cell +0.074 → +0.048 → +0.031 → +0.014), coverage 0.46 → 0.90 — **convergent, no bias floor**; the one-step early-t bias is purely an initialisation artifact. Production-scale arm: `mi_iter3` × 500 reps on all 24 core cells (sparse = headline, dense = no-harm check) + `mi_iter5` × 500 reps on the 8 rate-0.3 cells (convergence exhibit), methods indices 8/9 so all existing results stay valid. Spokes/CR out of scope for the extension (`pamm_ic_cr` has no `iter` yet). |

| 2026-06-11 | **Extension gate (council): conditional sign-off, fixes applied.** Codex + Claude legs verified all Fig-7 numbers against raw; Gemini signed off. Chained per-chain re-imputation confirmed methodologically sound (PMDA/Wei-Tanner/MICE-maxit family; chain refits do not invalidate Rubin pooling; surv se_ratio 0.99-1.07 at iter 5). Required fixes, all applied: (1) hazard scoping corrected — the sparse peaked-cell hazard gap at iter 5 (0.829-0.926) lies BELOW the oracle (0.911-0.925 in the same cells), i.e. it is IC information loss, unlike the dense-cell smoothing ceiling; (2) new disclosure — iteration occasionally amplifies weakly identified chains in flexible-TV n=300 models (extreme logHR reps, one degenerate hazard fit, zero mgcv warnings; iter-5 logHR RMSE 2.2-3.4x one-step even after excluding blow-ups; coverage protected only by exploding widths) — added to Fig 7, conclusions, and ?pamm_ic; (3) "cure"/"full nominal" wording weakened to computed levels (iter-5 surv min 0.927; rate-0.6 near-nominal within ±2pp); (4) aggregate-iter.R gains the Gate-R3 robustness kit (n_extreme blow-up guard excluding |est|>1e6 from moment summaries, pointwise-then-median se_ratio). Also fixed: error-row method labels in run-iter-pilot.R; "no longer depends on initialiser" softened to "progressively attenuated" (finite iter). |

ADEMP is locked at Gate R2 (pre-production); any later change goes in this log.
