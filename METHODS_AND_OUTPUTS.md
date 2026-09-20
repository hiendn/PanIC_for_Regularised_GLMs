# Locked methods and outputs

## Inherited simulation design

The study retains seven data-generating settings, a 20-dimensional coefficient
vector with ten active slopes, the inherited signal scaling, an independent
test sample of size 2,000, and a radius interval `[0,20]`. The primary path has
121 equally spaced radii. Calibration uses 31 log-spaced multipliers from
`0.01` to `100` and five balanced half-splits scored in both directions. The
solver controls and failure policy are shared across methods.

The frozen sources in `reference_implementation/` document the preceding
numerical engine. The active method set, support rule, seeds, and outputs are
defined by the current top-level configuration and simulation files.

## Assessed methods

Exactly three methods are assessed in every replication.

### PanIC-CF

For each calibration direction, the training-half GLM supplies coefficient
signs and the independent validation-half GLM supplies coefficient magnitudes:

```text
C_raw = sum_j sign(beta_training[j]) * beta_validation[j].
```

`C_raw` is retained for diagnostics and projected onto `[0,20]` for scoring.
Let `psi_hat` denote a candidate selected radius normalized to `[0,1]`, and let
`psi_target` be the normalized projected target. The calibration discrepancy is

```text
[psi_hat - psi_target]_+^2
  + sqrt(log(log(n_validation + exp(exp(1)))))
      * [psi_target - psi_hat]_+^2.
```

For each multiplier, the discrepancy is averaged over the ten directions. The
selected multiplier is the largest grid value whose mean is within one
descriptive row-level standard error of the minimum. This is an algorithmic
stabilization rule, not a confidence interval. The selected multiplier is then
used with the full-sample path.

### CV

An independently seeded balanced five-fold partition is used. At each radius,
the validation loss is averaged over folds. CV selects the smallest-indexed
radius attaining the minimum mean loss:

```text
j_CV = min argmin_j mean_fold loss(fold, j).
```

### BIC-like

The BIC-like comparator uses the full-sample path and a monotone cumulative
active count. The Gaussian criterion uses the inherited `log(n)/n` scaling;
the logistic and Poisson rows use the inherited one-half scaling and remain
exploratory analogues outside the Gaussian proposition. A small radius-index
term breaks criterion ties deterministically.

## Exact support rule

For every method,

```text
selected_j = (beta_hat_j != 0).
```

False positives, false negatives, exact recovery, and total support error use
this literal fitted-nonzero rule. The BIC-like active count is likewise
`sum_j(beta_hat_j != 0)` at each path point before the cumulative maximum is
taken. No numerical coefficient threshold is applied. The separate path-radius
interpolation tolerance remains a solver-accuracy diagnostic and is not a
support threshold.

## Primary estimands and joint decision

For scenario `s` and replication `r`, the paired support contrast is

```text
d_sr = (FP + FN)_PanIC-CF,sr - (FP + FN)_CV,sr.
```

If `dbar_s` and `se_s` are its scenario mean and Monte Carlo standard error,
the equal-weight effect and MCSE are

```text
Delta_support = (1/7) * sum_s dbar_s,
MCSE_support  = (1/7) * sqrt(sum_s se_s^2).
```

Support superiority passes when

```text
Delta_support + qnorm(0.95) * MCSE_support < 0.
```

The paired prediction contrast is

```text
q_sr = (D_PanIC-CF,sr - D_CV,sr) / D_CV,sr.
```

Every CV denominator must be finite and strictly positive. A violation stops
analysis and verification; no row is deleted and no denominator is floored or
replaced. Prediction noninferiority passes when the one-sided 95% upper Monte
Carlo bound for the equal-weight mean is below `0.001`. Margins `0.0005`,
`0.0025`, and `0.005` are sensitivity summaries only.

The joint conclusion requires complete PanIC-CF/CV pairing in all seven
settings and both component decisions to pass. Scenario-specific contrasts,
comparisons with BIC-like, and the grid study cannot rescue a failed joint
decision.

## Failure and numerical policy

If the full-sample path fails, all three method rows are marked failed. If
calibration fails while the full path remains valid, PanIC-CF uses the
prespecified multiplier `kappa=1` and records the failure and default. If a CV
fold path fails, the CV row is marked failed. No failed estimate is silently
replaced.

The production decision requires all 1,000 paired observations in every
scenario. The solver, interpolation, warning, projection, endpoint, and timing
diagnostics are retained even when they are not displayed in the compact
manuscript tables.

## Seed families

The fresh locked seeds are:

- main production: `2136092001`;
- grid production: `2138092001`;
- main smoke: `2140092001`; and
- grid smoke: `2142092001`.

The new method set and exact-nonzero policy were locked before these seed
families were used. Earlier development, production, grid, and smoke families
are listed as prior or retired and are included in the non-overlap audit.

## Main outputs

- `*_primary.csv`: three method rows per data-set replication.
- `*_diagnostics.csv`: path, calibration, CV, warning, endpoint, projection,
  multiplier, and timing diagnostics.
- `*_calibration_rows.csv`: ten calibration-direction rows with targets and
  the active PanIC-CF weight.
- `*_raw.rds`: lossless per-replication objects.
- `seed_ledger.csv`: every explicit training, calibration-split, CV-fold, and
  test seed.
- `simulation_summary.csv`: means and Monte Carlo standard errors for the
  three methods in seven scenarios.
- `paired_method_contrasts.csv`: the three pairwise method comparisons in each
  scenario.
- `paired_panic_cf_vs_cv.csv`: the PanIC-CF/CV rows used by the joint decision.
- `scenario_primary_estimands.csv` and `confirmatory_decision.csv`: the seven
  scenario effects and equal-weight joint gate.
- `prediction_margin_sensitivity.csv`: the primary and sensitivity margins.
- `relative_deviance_denominator_audit.csv`: the denominator-domain and
  finite-moment audit.
- `diagnostic_summary.csv` and `calibration_target_summary.csv`: compact
  numerical and target diagnostics.

## Grid outputs

The grid study uses 500 common-random-number replications at 61, 121, and 241
radii. Each grid has the same three assessed methods. The runner checks that
data, random streams, and raw calibration targets agree across grids before it
finishes.

The retained grid outputs comprise the replication files,
`grid_sensitivity_summary.csv`, paired across-grid contrasts, within-grid
method contrasts, the locked grid configuration, and
`table_grid_sensitivity.tex`. The summary has nine grid-by-method rows. Runtime
is not a column of the manuscript grid table.

## Deterministic manuscript tables

`Manuscript_Table_Rendering.R` renders the seven table mirrors from compact CSV
summaries without fitting a model or drawing a random number. The primary
performance table retains selected radius and signed attained-radius error but
does not repeat attained radius. The calibration table displays target, bias,
and multiplier summaries; detailed operational counts remain in the diagnostic
files. The grid table contains only the three active methods and omits runtime.

`manuscript/main.tex` embeds exact copies of the rendered tables between stable
marker comments. `Manuscript_Table_Tools.R --verify` checks byte-for-byte
agreement with the mirrors under `results/` and rejects external table-file
references.
