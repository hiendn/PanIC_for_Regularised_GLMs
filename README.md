# Consistent information criteria for regularised generalised linear models

R code, verification checks and manuscript-facing simulation outputs for the paper by Qingyuan Zhang and Hien Duy Nguyen.

## Requirements

The reference results use R 4.4.0, `glmnet` 4.1-10 and `Matrix` 1.7-0 on `aarch64-apple-darwin20`. Install the direct external dependency if needed:

```r
install.packages("glmnet")
```

## Reproduce and verify

Run the quick deterministic checks from the repository root in a fresh checkout:

```sh
Rscript tests/implementation_regression.R
Rscript Repository_Verification.R
```

The first command checks the split construction, seed separation, cross-signed target, projection, calibration weight, one-standard-error rule, CV allocation and deterministic completion of a full replication. The second verifies the schemas and numerical claims in the compact committed outputs. Both stop if a check fails.

Run the complete seven-setting confirmatory study and its replication-level checks with:

```sh
Rscript Simulation_Study.R --n-rep=500 --cores=4 --output-dir=results --overwrite
```

The study uses five independently seeded balanced half-splits in both directions, giving ten PanIC-CF calibration rows per data set. Each row combines signs from an ordinary training-half GLM with coefficient estimates from an independent ordinary validation-half GLM. The raw target is retained, projection onto `[0,20]` is used only for scoring, and calibration failure invokes the prespecified multiplier `kappa=1` with an explicit flag. The comparator is independently partitioned ordinary five-fold cross-validation.

The confirmatory master seed is `2026091802`. For setting index `s` and replication `r`, the training-data seed is `2026091802 + 100000*s + r`; the five split seeds add `10000 + 1000*(q-1)`, the CV seed adds `30000`, and the independent-test seed adds `60000`. Each worker receives these explicit seeds and does not generate its own random stream.

Run the supplementary numerical studies with:

```sh
Rscript Grid_Sensitivity_Study.R --n-rep=500 --cores=4 --output-dir=results
Rscript Runtime_Benchmark.R --output-dir=results
Rscript Boundary_Subsampling_Study.R
Rscript Solver_Validation.R
Rscript Supplemental_Results_Verification.R --output-dir=results
```

The grid study uses the disjoint master seed `2036091802` and common random numbers at 61, 121 and 241 radii. The boundary study uses 40,000 replications at each displayed sample size and 10,000 repeated subsamples for each subsample size. Runtime values are machine dependent and should not be expected to reproduce byte-for-byte.

Generated replication-level files, seed ledgers, timing records and checkpoints remain under `results/` but are ignored by Git. Use a fresh checkout or move existing generated files before changing the code or configuration.

## Source files

| File | Purpose |
| --- | --- |
| `Simulation_Study.R` | Main confirmatory entry point |
| `Simulation_Config.R` | Locked scenarios, grids, seeds and solver settings |
| `Regularised_GLM_Path_Functions.R` | Data generation, GLM loss evaluation and radius-path fitting |
| `PanIC_CF_Functions.R` | Repeated two-way PanIC-CF calibration and comparator logic |
| `Run_Confirmatory_Simulations.R` | Parallel confirmatory replications and seed ledger |
| `Simulation_Results_Analysis.R` | Monte Carlo summaries and manuscript tables |
| `Grid_Sensitivity_Study.R` | Common-random-number radius-grid comparison |
| `Runtime_Benchmark.R` | One-worker timing study |
| `Render_Runtime_Figure.R` | Runtime figure generation |
| `Implementation_Validation.R` | Deterministic implementation checks used by the main runner |
| `Results_Verification.R` | Replication-level confirmatory-result checks |
| `Supplemental_Results_Verification.R` | Grid and runtime checks |
| `Repository_Verification.R` | Checks of the compact committed outputs |
| `Auxiliary_Config.R` | Locked boundary and solver-check configuration |
| `Auxiliary_Simulation_Functions.R` | Auxiliary data-generation and constrained-path routines |
| `Boundary_Subsampling_Study.R` | Scalar boundary and subsampling study |
| `Solver_Validation.R` | Direct-fit and fixed-shape Gamma diagnostics |
| `tests/implementation_regression.R` | Standalone deterministic regression test |

## Reference outputs

The numerical implementation used for the retained reference outputs is identified by source commit [`7ab3f9b`](https://github.com/hiendn/PanIC_for_Regularised_GLMs/commit/7ab3f9b44486c8abbb45d3f847c8acbcbd2167ce).

The reference run contains 3,500 confirmatory data sets and 35,000 calibration rows. It recorded no full-path, calibration or CV failure; no default use; no solver warning; no target projection; and no multiplier or selected-radius endpoint selection. All 14 deterministic implementation checks in `results/verification_checks.csv` passed.

Relative to CV, the paired mean PanIC-CF total-support-error differences were `-0.472`, `-0.428`, `0.190`, `0.090`, `-0.426`, `0.050` and `-0.360` in the seven listed settings. PanIC-CF therefore reduced this support summary in the Gaussian and Poisson settings and had slightly higher error in the logistic settings. Every paired mean test-deviance difference had absolute value at most `1.13e-4` and was smaller than its own Monte Carlo standard error. This does not establish predictive equivalence; no systematic predictive difference was resolved at the achieved Monte Carlo precision.

The Gaussian active-count criterion is labelled BIC-like. Its logistic and Poisson counterparts are exploratory. The fixed-shape Gamma diagnostic emitted one nonconvergence warning, so Gamma is not included in the comparative simulation results.

| Manuscript item | Retained files |
| --- | --- |
| Primary support and performance tables | `results/simulation_summary.csv`, `results/paired_panic_cf_vs_cv.csv`, `results/table_primary_support.tex`, `results/table_primary_performance.tex` |
| Calibration table | `results/calibration_target_summary.csv`, `results/diagnostic_summary.csv`, `results/table_calibration.tex` |
| Grid-sensitivity table | `results/grid_sensitivity_summary.csv`, `results/grid_sensitivity_paired_contrasts.csv`, `results/table_grid_sensitivity.tex` |
| Boundary and subsampling table/figure | `results/boundary_summary.csv`, `results/subsampling_summary.csv`, `results/table_boundary.tex`, `results/table_subsampling.tex`, `results/figure_boundary_subsampling.pdf` |
| Runtime figure | `results/runtime_summary.csv`, `results/figure_runtime.pdf` |
| Solver and Gamma diagnostics | `results/solver_validation.csv`, `results/gamma_validation.csv` |

## Contact

Hien Duy Nguyen: h.nguyen5@latrobe.edu.au
