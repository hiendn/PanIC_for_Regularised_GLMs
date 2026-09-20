# Consistent information criteria for regularised generalised linear models

R code, verification checks, and manuscript-facing simulation outputs for the
paper by Qingyuan Zhang and Hien Duy Nguyen.

The active confirmatory implementation assesses exactly three methods:

- `PanIC-CF`, with the square-root log-log calibration weight;
- `CV`, ordinary five-fold cross-validation selecting the smallest radius that
  attains the minimum mean validation loss; and
- `BIC-like`, the family-appropriate active-count comparator.

The vendored files under `reference_implementation/` preserve the numerical
lineage of the project. They supply shared historical primitives but do not
define an additional assessed method or the active support rule.

## Confirmatory scope

The main study has seven equally weighted Gaussian, logistic, and Poisson
settings and 1,000 paired replications per setting. Each replication produces
three method rows. The joint comparison of PanIC-CF with CV requires both:

1. support-error superiority, assessed by the one-sided 95% upper Monte Carlo
   bound for the paired difference in `FP + FN`; and
2. predictive noninferiority, assessed by the corresponding upper bound for
   `(D_PanIC-CF - D_CV) / D_CV` against the fixed margin `0.001`.

All seven scenarios must be completely paired. Scenario-specific and grid
comparisons are descriptive and cannot replace either joint endpoint.

For the committed production run, the equal-weight PanIC-CF-minus-CV support
contrast was `-0.421857` (MCSE `0.018212`; one-sided 95% upper bound
`-0.391901`). The corresponding relative test-deviance contrast was
`0.000119699` (MCSE `0.000032143`; upper bound `0.000172569`). The first upper
bound is below zero and the second is below the prespecified `0.001` margin, so
the bounded joint conclusion passed for these seven simulation settings. This
is a finite-design result, not a claim that PanIC-CF dominates CV generally.

Support is the literal nonzero pattern of the fitted coefficient vector:
`beta_hat != 0`. The same exact-nonzero rule supplies the active count in the
BIC-like criterion. No post hoc coefficient threshold is used.

The method set and exact-zero rule were locked before the fresh seed families
were used for regeneration. This statement concerns the new regeneration; it
does not imply that these revisions preceded earlier exploratory or historical
outputs.

## Requirements

The preceding production environment used R 4.4.0, `glmnet` 4.1-10, and
`Matrix` 1.7-0 on `aarch64-apple-darwin20`. The only direct external package
dependency is `glmnet`:

```r
install.packages("glmnet")
```

The scripts otherwise use base R. A repository-local `Rlib/` may be used and
is ignored by Git.

## Validation

From the repository root, run:

```sh
Rscript Implementation_Validation.R
Rscript Repository_Verification.R
```

The deterministic implementation checks cover the method lock, exact-nonzero
support and BIC counts, balanced splitting, ordinary minimum-loss CV,
end-to-end fitting, and separation of production, smoke, and retired seed
streams. The repository verifier checks the committed schemas and decisions,
reconstructs the displayed summaries, and compares the rendered table mirrors
with the embedded manuscript tables. The main and grid pipelines also invoke
their replication-level verifiers.

## Smoke reproduction

Smoke mode uses seed families reserved for diagnostics and cannot write to the
committed `results/` directory:

```sh
Rscript Simulation_Study.R --smoke --n-rep=2 --cores=4 \
  --output-dir=smoke_results --overwrite
Rscript Grid_Sensitivity_Pipeline.R --smoke --n-rep=2 --cores=4 \
  --output-dir=smoke_results_grid
```

The dedicated smoke seeds are `2140092001` for the main pipeline and
`2142092001` for the grid pipeline.

## Full reproduction

The main regeneration uses master seed `2136092001`, seven settings, and 1,000
replications per setting:

```sh
Rscript Simulation_Study.R --n-rep=1000 --cores=4 \
  --output-dir=production_results --overwrite
```

The grid regeneration uses the disjoint master seed `2138092001`, 500
common-random-number replications, and 61, 121, and 241 equally spaced radii on
`[0,20]`:

```sh
Rscript Grid_Sensitivity_Pipeline.R --n-rep=500 --cores=4 \
  --output-dir=production_results
```

The runners record the exact configuration, explicit seed ledger, environment,
replication-level rows, compact summaries, validation records, and manuscript
table mirrors. Release manifests and checksums are refreshed only after the
new runs and all verifiers complete.

## Outputs and table presentation

The main raw files contain three method rows per data-set replication, one
diagnostic row, ten calibration-direction rows, and a lossless R object. The
principal compact outputs are:

- `simulation_summary.csv` and `paired_method_contrasts.csv`;
- `paired_panic_cf_vs_cv.csv`, `scenario_primary_estimands.csv`, and
  `confirmatory_decision.csv`;
- `relative_deviance_denominator_audit.csv` and
  `prediction_margin_sensitivity.csv`;
- `diagnostic_summary.csv` and `calibration_target_summary.csv`; and
- the corresponding grid summaries and paired contrasts.

The displayed tables are intentionally compact. The primary performance table
reports the signed attained-radius error and selected grid radius, rather than
also repeating the attained radius. The calibration table keeps the target,
bias, and multiplier summaries; sparse operational counts remain in the
diagnostic outputs and are stated in prose when relevant. The grid table shows
the three assessed methods at each grid size and omits runtime. Runtime is
reported separately using the method labels `PanIC-CF`, `CV`, and `BIC-like`.

`manuscript/main.tex` embeds all seven numerical table environments directly
between stable marker comments. The files under `results/table_*.tex` are
deterministic, checksum-controlled mirrors rather than manuscript inputs. They
can be regenerated without fitting a model:

```sh
Rscript Render_Manuscript_Tables.R \
  --results-dir=results --output-dir=/tmp/panic-tables
Rscript Manuscript_Table_Tools.R --verify \
  --manuscript=manuscript/main.tex --tables-dir=results
```

## Reproducibility safeguards

- The method, estimands, margins, failure policy, and seed families are
  recorded in `METHOD_LOCK_SHA256.txt`.
- The production seeds are `2136092001` and `2138092001`; the dedicated smoke
  seeds are `2140092001` and `2142092001`.
- Prior, development, and retired production/smoke seed families are excluded
  explicitly and checked for stream overlap.
- Every training, calibration-split, CV-fold, and test seed is a deterministic
  function of its master seed, scenario, and replication.
- A calibration failure uses the prespecified `kappa=1` default and records the
  failure and default flag.
- Every CV denominator in the relative-deviance endpoint must be finite and
  strictly positive. A violation stops analysis; no replication is removed and
  no denominator is replaced.
- Logistic and Poisson BIC-like rows remain exploratory active-count analogues
  outside the Gaussian proposition.

See `METHODS_AND_OUTPUTS.md` for the formal estimands and output schema and
`PROVENANCE.md` for the revision, seed, and packaging history.

## Main source files

| File | Purpose |
| --- | --- |
| `Simulation_Config.R` | Locked scenarios, grids, seeds, margins, and numerical settings |
| `PanIC_CF_Functions.R` | PanIC-CF calibration and the three-method replication driver |
| `Run_Confirmatory_Simulations.R` | Seven-setting runner and seed ledger |
| `Simulation_Results_Analysis.R` | Summaries, paired endpoints, decision, and audits |
| `Results_Verification.R` | Independent main-result reconstruction |
| `Grid_Sensitivity_Study.R` | Fresh-seed common-random-number grid study |
| `Grid_Results_Verification.R` | Independent grid-result reconstruction |
| `Implementation_Validation.R` | Method, support-rule, and seed checks |
| `Repository_Verification.R` | Compact-release and table verification |
| `Runtime_Benchmark.R` | Three-method one-worker timing benchmark |
| `Manuscript_Table_Rendering.R` | Deterministic renderers for the seven table mirrors |
| `reference_implementation/` | Frozen historical numerical lineage, not an assessed method set |
| `manuscript/` | Self-contained manuscript and reviewer response |

## Citation and contact

Citation metadata are in `CITATION.cff`; no DOI has been assigned in this
repository.

Hien Duy Nguyen: h.nguyen5@latrobe.edu.au
