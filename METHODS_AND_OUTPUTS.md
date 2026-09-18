# Locked methods and outputs

## Inherited design

The seven data-generating settings, 20-dimensional coefficient vector with ten
active slopes, signal scaling, independent test size 2,000, radius interval
`[0,20]`, 121-point primary radius grid, 31-point multiplier grid, five
balanced half-splits scored in both directions, solver tolerances, and failure
policy are inherited without change from the preceding locked production
release.

For a calibration direction, the training-half GLM supplies coefficient signs
and the independent validation-half GLM supplies coefficient magnitudes. The
raw cross-signed target is

```text
C_raw = sum_j sign(beta_training[j]) * beta_validation[j]).
```

The raw value is retained for diagnostics and projected onto `[0,20]` only for
scoring. The selected multiplier is the largest member of the unchanged grid
whose mean discrepancy is within one descriptive row-level standard error of
the minimum over the ten directions. This is an algorithmic stabilization
rule, not a confidence interval.

## Revised and original PanIC-CF weights

Let `psi_hat` be the normalized radius selected for a calibration path and
`psi_target` the projected normalized target. Revised PanIC-CF uses

```text
[psi_hat - psi_target]_+^2
  + sqrt(log(log(n_V + exp(exp(1)))))
      * [psi_target - psi_hat]_+^2.
```

The same-path sensitivity `PanIC-CF-original` replaces the square-root
log-log weight by the preceding log-log weight. Both discrepancies are formed
during the same loop over the same selected calibration radii. Each has its
own one-standard-error multiplier selection, but the sensitivity requires no
additional path, pilot, or test fit.

The slower primary weight still diverges and is `o(n)`, so the existing
asymmetric-target consistency condition is retained. The multiplier remains
in the same fixed positive compact grid, so the uniform-multiplier PanIC
argument is unchanged.

## CV comparators

The independently seeded balanced five-fold partition is unchanged. At every
radius, the mean and standard error of the five fold-specific validation means
are calculated. `CV-min` chooses the smallest-indexed minimizer of mean loss,
which exactly reproduces the preceding comparator.

`CV-1SE` chooses

```text
min {j : mean_loss[j] <= mean_loss[j_min] + SE_loss[j_min]},
```

where radii are in strictly increasing order and
`SE_loss[j] = sd(fold_loss[,j])/sqrt(5)`. Thus it chooses the smallest, most
regularized eligible radius. It reuses the CV-min fold paths and full-sample
path and is strictly secondary.

## Primary estimands and joint decision

For scenario `s`, let `d_sr` be the paired total-support-error contrast in
replication `r`:

```text
d_sr = (FP + FN)_PanIC-CF,sr - (FP + FN)_CV-min,sr.
```

If `dbar_s` and `se_s` are its scenario mean and Monte Carlo standard error,
the locked equal-weight primary effect and MCSE are

```text
Delta_support = (1/7) * sum_s dbar_s,
MCSE_support  = (1/7) * sqrt(sum_s se_s^2).
```

The scenario streams are disjoint. Support superiority passes when

```text
Delta_support + qnorm(0.95) * MCSE_support < 0.
```

For prediction, the replication-level paired contrast is scale-free:

```text
q_sr = (D_PanIC-CF,sr - D_CV-min,sr) / D_CV-min,sr.
```

The equal-weight mean and MCSE use the same formulas. Prediction
noninferiority passes when its one-sided 95% upper normal Monte Carlo bound is
below the prospectively fixed margin `0.001`. Sensitivity margins `0.0005`,
`0.0025`, and `0.005` are reported without changing the decision.

Complete PanIC-CF/CV-min pairing in all seven settings is also required. The
support and prediction requirements form an intersection-union decision: a
joint claim is made only if both one-sided level-0.05 component requirements
pass. No multiplicity reduction is needed to control the union null when both
components are required.

Scenario-specific contrasts, comparisons involving PanIC-CF-original or
CV-1SE, and the grid study are secondary. They cannot rescue a failed locked
decision and do not justify universal dominance over cross-validation.

## Failure and numerical policy

If a full-sample path fails, every method depending on that path is marked
failed. If calibration alone fails while the full path is valid, the unchanged
prespecified multiplier `kappa=1` is used and both calibration failure and
default use are recorded. If a CV fold path fails, both CV-min and CV-1SE are
marked failed. No failed estimate is silently replaced.

The joint claim requires all 1,000 paired observations in all seven settings;
otherwise `complete_pairing=0` and `joint_claim_pass=0`, even if the numerical
bounds happen to pass.

## Main outputs

- `*_primary.csv`: five method rows per data-set replication.
- `*_diagnostics.csv`: failures, warnings, multipliers, CV indices and
  one-standard-error thresholds, projection counts, and interpolation error.
- `*_calibration_rows.csv`: ten direction rows with revised and original
  weights.
- `*_raw.rds`: lossless per-replication R objects.
- `seed_ledger.csv`: every explicit training, split, CV, and test seed.
- `simulation_summary.csv`: method means and Monte Carlo standard errors.
- `paired_method_contrasts.csv`: six prespecified paired method comparisons.
- `scenario_primary_estimands.csv`: seven support and relative-deviance
  effects with one-sided upper bounds.
- `confirmatory_decision.csv`: the locked aggregate effects and joint gate.
- `prediction_margin_sensitivity.csv`: primary and sensitivity margins.
- `diagnostic_summary.csv` and `calibration_target_summary.csv`: numerical and
  target diagnostics.
- `manuscript_generated/table_primary_support.tex` and
  `table_confirmatory_decision.tex`: manuscript-ready tables.

## Grid outputs

The fresh-seed grid study produces replication files for 61, 121, and 241
radii, `grid_sensitivity_summary.csv`, paired across-grid contrasts,
within-grid method contrasts, a locked grid configuration, and
`manuscript_generated/table_grid_sensitivity.tex`. Data, streams, raw targets,
and both weights must match across grids before the runner completes.
