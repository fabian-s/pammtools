# Interval-Censored Survival Support for `pammtools` via Multiple Imputation

## Context

`pammtools` currently transforms only right-censored and left-truncated data into
PED format (the `as_ped` → `split_data` → `survival::survSplit` path supports 2-arg
and 3-arg `Surv` only). Interval-censored (IC) data — where the event time for
subject *i* is known only to lie in `(L_i, R_i]` — has **no support**. The JSS paper's
Limitations section (and Fabian's inline TODO) flags this, noting that folding IC
*directly* into the PEM/Poisson likelihood while preserving Poisson proportionality
is hard/maybe impossible.

The defensible way around this is **multiple imputation (MI) + refit** (cf. Pan 2000
for Cox, Delord & Génin 2016 for competing risks): an *imputed* event time is an
*exact* event time, so once we draw `T_i ∈ (L_i, R_i]` from the model-based
conditional distribution, the **entire existing right-censored pipeline is reused
unchanged**. We repeat this `M` times and **pool the posterior draws** through the
package's existing simulation-based inference machinery. This delivers calibrated
inference (naive midpoint imputation invalidates CIs even when point estimates are
okay) with very little new numerical code.

**Scope (after narrowing with the user):**
- Implement **single-event AND competing risks** IC support now; recurrent/MSM =
  assessed future work only.
- Default imputation = **model-based conditional hazard** draw (midpoint only as
  initializer).
- Pooling = **pooled posterior draws** (reuse existing `get_sim_ci*`/`get_cif`).
- **Dropped from scope** (judged low value-for-effort): Rubin's-rules closed-form
  pooling and per-functional within-variance plumbing; the analytic Louis-style
  missing-information / variance-inflation approach. *Note:* a fraction-of-missing-
  information (FMI) diagnostic is cheaply recoverable later from the pooled draws
  (between-fit vs within-fit variance) if ever wanted — no new infrastructure needed.

**Intended outcome:** users specify `Surv(L, R, type = "interval2")`, call a
`pamm_ic()` driver, and get the full `add_hazard` / `add_surv_prob` / `add_cif` /
RMST family with MI-valid confidence intervals, plus a simulation harness for
coverage studies.

---

## Cross-cutting invariant (the one correctness trap)

**`cut` is resolved exactly once and shared across all `M` imputations.** If breaks
were re-derived per imputation (the `get_cut(cut = NULL)` default keys off observed
event times), the interval structure would differ across fits and pooling would be
ill-defined. The IC pipeline freezes `cut` in `as_ped_ic`/`pamm_ic` and passes it to
every `split_data`/`as_ped_cr` call.

> **Implementation correction (verified empirically, R available).** The plan
> originally assumed a fixed `cut` also makes the design matrix
> `predict.gam(fit_m, newdata, "lpmatrix")` *identical* across fits, so pooled draws
> could share one `X`. **This is false:** mgcv's identifiability centering (and
> thin-plate knot placement) makes the *constrained* spline basis depend on each
> imputed data set, so `X` differs across fits (`max|X_1 - X_m| > 0`, confirmed).
> The pooling layer was therefore built to propagate **each fit's draws through its
> own design matrix** (`ic_ci_draws`, `ic_cif_draws_group` in `R/pool-ic.R`) and pool
> the resulting functional draws — the statistically correct approach, needing no
> shared-`X` assumption. Point estimates are the average of per-fit point predictions
> (the MI estimate). The sampler's calibration is verified by a probability-integral-
> transform KS test; a coverage study confirms MI widens intervals and restores
> coverage relative to naive midpoint imputation.

---

## Implementation steps

### 1. Input representation & preprocessing — `R/as-ped-ic.R` *(risk: low)*
- `detect_ic(formula, data)`: evaluate the LHS `Surv(...)` in `data`, inspect
  `attr(., "type")`. Return `"interval2"` when type ∈ {`interval`,`interval2`};
  `counting` (left-trunc) and 2-arg right-censoring fall through untouched →
  **zero behavior change for existing users**.
- Intercept in `as_ped.data.frame` (`R/as-ped.R:151`) **before** `status_error()`
  (line 163, which would choke on survival's 0/1/2/3 interval status codes) and the
  existing `transition`/`get_event_types` dispatch. Route IC to `as_ped_ic()`.
- `parse_ic_surv(formula, data)` → tibble `(ic_L, ic_R, ic_kind, <covars>)` with
  `ic_kind ∈ {exact, right, left, interval}` from survival's status codes
  (`1`=exact `L==R`, `0`=right `R=Inf`, `2`=left `L=0`, `3`=interval).
- `resolve_ic_cut(ic, cut, max_time)`: if `cut` supplied, freeze it; else default to
  the sorted unique **finite interval endpoints** (∪ of `L`,`R`) capped at `max_time`
  — statistically natural (breaks at inspection times, Turnbull-like).
- `as_ped_ic()` returns the **midpoint-initializer PED** (class
  `c("ped_ic_init","ped",...)`) with the parsed bounds and frozen `cut` carried as
  `attr(ped, "ic")` / `attr(ped, "breaks")`, so `as_ped` always yields a usable `ped`.
- Per-imputation PED build is just the existing 2-arg path:
  `split_data(Surv(.t_imp, .ev) ~ covars, imp_data, cut = <frozen>)`. Reuses
  `split_data` (`R/split-data.R:14`) verbatim. Edge cases to respect: clamp `t_imp`
  strictly inside `(L,R]` and `< max(cut)` to avoid the `tstart == ped_time` filter
  (`split-data.R:136/148`) and the `ped_status` zeroing at `max(cut)` (lines 124/141).

### 2. Conditional-hazard sampler — `R/impute-ic.R` *(risk: low–med)*
- `impute_ic_times(object, ic, cut, beta = NULL, newdata_template = NULL)` → numeric
  `t_imp`. Exact/right rows pass through unchanged (right-censored subjects are **not**
  imputed — they contribute correctly as censored at `L`).
- **Math (inverse-CDF on the piecewise-linear cumulative hazard):** with `H_i(t)`
  piecewise-linear (rates `h_{ij}` on the fixed grid), draw `U~Unif(0,1)` and set
  `Λ* = H(L) − log(1 − U·(1 − exp(−(H(R) − H(L)))))`, then invert
  `t = c_{j*−1} + (Λ* − H(c_{j*−1})) / h_{i,j*}` where `j*` brackets `Λ*`.
  - Left-censored: `H(L)=0`. Interval: general form. Degenerate `H(R)−H(L) < eps`:
    Taylor branch → uniform draw `L + U·(R−L)`.
- **Vectorized hazard evaluation (fixed-grid accumulation, reuses package machinery):**
  build one per-`(subject, interval)` grid on the frozen `cut` via `make_newdata` /
  `reconstruct_intlen` / `reconstruct_cutpoints` (`R/make-newdata.R`,
  `R/add-functions.R`); one `predict.gam(object, grid, "lpmatrix")` → cache `X`;
  `h = exp(X %*% beta)`; grouped `cumsum` for `H_grid`; `findInterval` to locate
  `L,R,Λ*`. `s(tend)`, `te(tend,x)`, `by=` terms are handled automatically since they
  key off `tend`. Cache `X` once and only recompute `h` per imputation.

### 3. MI driver loop — `R/pamm-ic.R` *(risk: low–med)*
- `pamm_ic(formula, data, cut = NULL, max_time = NULL, m = 10L, proper = TRUE,
  init = c("midpoint","uniform"), engine = "gam", ...)`.
  1. Parse + freeze grid via `as_ped_ic`.
  2. Initializer: midpoint PED → `fit0 <- pamm(model_formula, ped0, engine, ...)`.
  3. Loop `mm = 1..m`: `beta_mm <- if (proper) MASS::mvrnorm(1, coef(anchor),
     anchor$Vp) else coef(anchor)` → `impute_ic_times(...)` → `split_data(...,
     cut = frozen)` → `pamm(...)`. Store fit.
- **Proper MI (`proper = TRUE`) is the default**: drawing `β^(m) ~ N(β̂, Vp)` before
  imputing is the Bayesianly-proper choice that makes pooled draws valid; `N(β̂, Vp)`
  is the *same* posterior `get_sim_ci` already uses, so it's internally consistent.
  `proper = FALSE` documented as diagnostics-only.
- Returns `structure(list(fits, init_fit, ic, cut, formula, m, proper),
  class = c("pamm_ic","list"))`. Each element of `fits` is a plain `pamm`, so every
  existing `add_*` works on any single fit. Add `print/summary.pamm_ic`.

### 4. Competing-risks extension — `R/pamm-ic.R` *(risk: med)*
- `pamm_ic_cr(formula, data, cut = NULL, m, proper, cause_unknown = NULL, ...)` and
  `impute_ic_cr(object, ic, cut, beta)` → `(t_imp, cause_imp)`.
- Imputation model = cause-specific PAMM on stacked `ped_cr` (`s(tend, by = cause)`).
  Draw `T` from the **all-cause** conditional hazard `H_•=Σ_k H_k` (step 2 sampler).
  - **Cause known** (timing IC, cause observed): keep cause `k`, draw `T` by
    rejection — propose from all-cause conditional, accept with prob
    `h_k(T)/h_•(T) ≤ 1`.
  - **Cause unknown**: after drawing `T`, sample cause with prob `h_k(T)/h_•(T)` —
    exactly the `sample(seq_rhs, prob = probs)` pattern at `sim-pexp.R:350` and the
    `hazard/total_hazard` ratio in `get_cif` (`add-functions.R:1010`).
- Build per-imputation PED via existing `as_ped_cr` (`R/as-ped.R:339`); **freeze the
  per-cause `cut` list** (it builds per-event-type cuts at line 353) and pass it
  explicitly each imputation. Reference: Delord & Génin (2016).

### 5. Pooling layer (pooled posterior draws) — `R/pool-ic.R` *(risk: low)*
- `pooled_sim_coef(pamm_ic, nsim = 500L)` stacks `rmvnorm(nsim/m, coef(f), f$Vp)`
  across fits → a sample from the MI mixture posterior `(1/M)Σ N(β̂^(m),Vp^(m))`,
  whose variance ≈ `W + B` (Rubin total var, modulo the `(1+1/M)` factor, addable by
  scaling the between-fit spread by `sqrt(1+1/M)`).
- Add `add_*.pamm_ic` S3 methods (`add_hazard`, `add_cumu_hazard`, `add_surv_prob`,
  `add_cif`, RMST) that build `X` from `fits[[1]]` (shared `X`!), get the pooled
  matrix, and call the **existing** propagation loops in `get_sim_ci` /
  `get_sim_ci_surv` / `get_cif` / trans-prob **verbatim**. Load-bearing inner loops
  untouched — this is the whole point: postprocessing is a thin adapter.
- Point estimate = `colMeans` of the stacked draws (equiv. average of per-fit
  predictions).
- *(Optional, deferred)* an FMI/between-vs-within diagnostic can be read off the same
  stacked draws later if wanted — no extra infrastructure.

### 6. Simulation & validation tooling — `R/sim-pexp.R` + `tests/` *(risk: low)*
- `add_inspections(sim_df, schedule = NULL, rate = 1, mechanism = c("random","fixed","mixed"))`:
  turn exact `sim_pexp` / `sim_pexp_cr` output into IC data by generating per-subject
  inspection times (e.g. cumulative `Exp(rate)`), setting `L` = last inspection
  `< time`, `R` = first inspection `≥ time` (`R=Inf` if beyond last). Keep true `time`
  as a hidden column for coverage scoring. CR variant records cause at the positive
  inspection (optionally hidden to test unknown-cause imputation).
- **Coverage study** (heavier script in `inst/`/vignette, not CRAN check): across
  scenarios (n, inspection density, baseline shape, effect size, censoring rate, K),
  compare naive-midpoint vs MI-pooled-draws on coverage, width, bias, RMSE; optionally
  cross-check `S(t)` against `icenReg` Turnbull NPMLE (`icenReg` a *Suggests* dep only).

### 7. Recurrent / MSM — roadmap only *(risk: high, deferred)*
Document in a vignette/section. The MI+refit skeleton extends (per-transition
sampler; PED via existing `split_data_multistate` / `as_ped_multistate`,
`R/split-data.R:213`), but two genuinely hard parts must be flagged: (1)
**sequential/joint imputation** — an imputed transition time is the left-truncation
entry for the next, so imputations are a path draw with ordering constraints, not
independent marginals; clock-forward vs clock-reset changes how this propagates; (2)
**congeniality** — independent per-transition imputation ignoring within-subject
dependence/frailty is uncongenial and biases variance.

---

## New artifacts summary

| Piece | New functions | File | Key reuse |
|---|---|---|---|
| Preproc | `detect_ic`, `parse_ic_surv`, `as_ped_ic`, `resolve_ic_cut` | `R/as-ped-ic.R` | `as_ped.data.frame` dispatch, `split_data` |
| Sampler | `impute_ic_times` | `R/impute-ic.R` | `predict.gam` lpmatrix, `make_newdata` |
| Driver | `pamm_ic`, `print/summary.pamm_ic` | `R/pamm-ic.R` | `pamm`, `split_data`, `Vp` |
| CR | `pamm_ic_cr`, `impute_ic_cr` | `R/pamm-ic.R` | `as_ped_cr`, `get_cif`/`sim_pexp_cr` cause logic |
| Pooling | `pooled_sim_coef`, `add_*.pamm_ic` | `R/pool-ic.R` | `get_sim_ci*`, `get_cif` |
| Sim/test | `add_inspections` + coverage harness | `R/sim-pexp.R`, `tests/` | `sim_pexp`, `sim_pexp_cr`, `rpexp` |

Follow package conventions: roxygen2 with `@export`/`@rdname`/`@keywords internal`,
`checkmate` assertions, `mvtnorm`/`MASS` for draws (already deps or trivial adds),
cached `tests/testthat/helper-fixtures.R`-style fixtures with `set.seed`.

---

## Verification

1. **Unit tests** (`tests/testthat/test-interval-censored.R`, follow existing style):
   - **Fixed-cut invariant**: `predict.gam(fit_m, nd, "lpmatrix")` identical across all
     `m` fits from a `pamm_ic` object (the load-bearing precondition).
   - Sampler: every `t_imp ∈ (L, R]`; left/right/exact/degenerate edge cases;
     empirical CDF of draws matches the analytic truncated PCH CDF (KS test).
   - Exact-only data (`L==R` for all) → `pamm_ic` reproduces a plain `pamm` fit.
   - CR: imputed causes' empirical split matches `h_k/h_•`; per-cause cut frozen.
   - Pooling: `S(t)` monotone; pooled total variance ≈ `W + (1+1/M)B`.
2. **Numerical/statistical validation** (script, not CRAN check):
   - `add_inspections` coverage study — confirm midpoint undercovers and
     MI-pooled-draws hits ~nominal 95%; cross-check `S(t)` against `icenReg` Turnbull.
3. **Regression**: full `devtools::test()` green (no behavior change on non-IC paths,
   guaranteed by `detect_ic` early-return); `R CMD check --as-cran` clean.
4. **Worked example**: end-to-end on IC-censored simulated (and/or `tumor`) data —
   `pamm_ic()` → `add_surv_prob()` / `add_cif()` — in a new
   `vignettes/interval-censored.Rmd`, paralleling `competing-risks.Rmd`.

## Suggested build order
1 → 2 → 3 (single-event MVP, testable end-to-end) → 6 (sim tooling, enables coverage)
→ 5 (pooled-draws postprocessing) → 4 (CR) → 7 (roadmap docs).
