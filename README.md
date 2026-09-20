# Consistent information criteria for regularised generalised linear models

R code, verification checks, and manuscript-facing simulation outputs for the
paper by Qingyuan Zhang and Hien Duy Nguyen.

This repository contains the self-contained reference implementation of
PanIC-CF used for the final confirmatory simulations. It also contains the
ordinary five-fold cross-validation comparators `CV-min` and `CV-1SE`, the
same-path `PanIC-CF-original` sensitivity analysis, the BIC-like active-count
comparator, and the retained boundary, subsampling, solver, and runtime
diagnostics.

## Confirmatory result and its scope

The primary comparison was fixed before the production results were inspected.
Across seven equally weighted Gaussian, logistic, and Poisson settings, with
1,000 paired replications per setting, PanIC-CF had lower total support error
than `CV-min`:

```text
PanIC-CF minus CV-min support error = -0.438143
Monte Carlo standard error         =  0.017999
one-sided 95% upper bound          = -0.408538
```

The paired relative test-deviance difference was `0.000159955` (about
`0.016%`), with Monte Carlo standard error `0.000032904` and one-sided 95%
upper bound `0.000214076`. This is below the prespecified noninferiority margin
`0.001` (0.1%). Complete pairing and both endpoint decisions were required, so
the locked joint decision passed.

The claim is deliberately narrower than “PanIC-CF dominates cross-validation.”
`CV-1SE` achieved substantially lower support error than PanIC-CF but about
4.9% worse test deviance in the equal-weight secondary comparison. Accordingly,
the confirmatory claim is support-error superiority to ordinary minimum-loss
five-fold CV, jointly with predictive noninferiority under the stated margin.
Scenario-specific, original-weight, `CV-1SE`, and grid comparisons are
secondary.

## Requirements

The production results used R 4.4.0, `glmnet` 4.1-10, and `Matrix` 1.7-0 on
`aarch64-apple-darwin20`. Install the direct external dependency if needed:

```r
install.packages("glmnet")
```

The scripts use base R and `glmnet`; no package manager or private source tree
is required. If desired, packages may be installed into a repository-local
`Rlib/` directory, which is ignored by Git.

## Verify the committed release

Run these checks from the repository root in a fresh checkout:

```sh
Rscript Implementation_Validation.R
Rscript Repository_Verification.R
```

`Implementation_Validation.R` performs 21 deterministic checks covering the
method lock, seed-stream separation, the revised and original weight formulas,
balanced splitting, CV-min/CV-1SE selection, deterministic end-to-end fitting,
and exact agreement of the inherited methods with the frozen preceding
implementation.

`Repository_Verification.R` checks the compact result manifest and schemas,
reconstructs the two equal-weight confirmatory endpoints from the seven
scenario rows, checks every displayed decision against the locked thresholds,
audits diagnostic counts, checks the grid study, regenerates all seven table
mirrors byte-for-byte, and verifies their exact embedding in the manuscript.
The production pipeline's independent `Results_Verification.R` reconstructs
every primary support/performance field used by the tables and the
relative-deviance denominator audit from replication-level rows.
`Grid_Results_Verification.R` likewise reconstructs every displayed grid
summary field from the raw grid rows before checking the rendered table.
Every validation script stops if a critical check fails.

## Fast smoke reproduction

Diagnostic runs must use the dedicated smoke streams. The following commands
run two replications in every main setting and two common-random-number
replications at each grid size, then invoke the full replication-level
verifiers:

```sh
Rscript Simulation_Study.R --smoke --n-rep=2 --cores=4 \
  --output-dir=smoke_results --overwrite
Rscript Grid_Sensitivity_Pipeline.R --smoke --n-rep=2 --cores=4 \
  --output-dir=smoke_results_grid
```

Smoke mode uses master seeds `2106091802` and `2126091802`. It cannot write to
the committed `results/` directory and cannot consume a production seed.

## Full reproduction

The confirmatory study uses main master seed `2066091802`, seven settings, and
1,000 replications per setting. Run it in a separate ignored directory:

```sh
Rscript Simulation_Study.R --n-rep=1000 --cores=4 \
  --output-dir=production_results --overwrite
```

The pipeline verifies the method lock, runs all seven settings, writes the seed
ledger and replication-level diagnostics, creates the Monte Carlo summaries
and manuscript tables, and independently reconstructs the joint decision from
the raw rows.

The grid study uses disjoint master seed `2086091802`, 500 common-random-number
replications, and 61, 121, and 241 equally spaced radii on `[0,20]`:

```sh
Rscript Grid_Sensitivity_Pipeline.R --n-rep=500 --cores=4 \
  --output-dir=production_results
```

The complete main and grid runs generate about 33 MB of replication-level
CSV/RDS output. Those regenerable files are ignored by Git. Their names, sizes,
and SHA-256 hashes at completion of the locked production run are recorded in
`results/production_file_manifest.csv`. That file is an immutable historical
production snapshot, so it retains the original `manuscript_generated` paths;
the self-contained release archive has its own payload manifest. The compact
summaries, configuration, seed ledger, canonical table mirrors, and
verification records are committed.

## Self-contained manuscript tables

`manuscript/main.tex` embeds all seven numerical table environments
directly between stable marker comments. It therefore has no `\input`,
`\include`, or `\IfFileExists` dependency on a `table_*.tex` file. The seven
files under `results/table_*.tex` are retained as deterministic,
checksum-controlled reproducibility mirrors of those embedded blocks, not as
manuscript compilation inputs. Their SHA-256 hashes are recorded in
`results/RESULTS_SHA256.txt`.

Regenerate the seven mirrors from the committed compact summaries without
rerunning a simulation:

```sh
Rscript Render_Manuscript_Tables.R \
  --results-dir=results --output-dir=/tmp/panic-tables
```

`Render_Manuscript_Tables.R` is the command-line wrapper for the pure renderers
in `Manuscript_Table_Rendering.R`. To compare the retained mirrors with the
committed manuscript source, run:

```sh
Rscript Manuscript_Table_Tools.R --verify \
  --manuscript=manuscript/main.tex --tables-dir=results
```

The optional `--refresh` mode replaces the seven marker-delimited blocks and
then performs the same byte-for-byte comparison. Verification also fails if
an external table-file reference remains. This table rendering and embedding
amendment changes presentation and reproducibility only; it does not rerun the
production study or change any reported numerical result.

## Reproducibility safeguards

- The statistical method, estimands, margins, failure policy, and production
  seeds are recorded in `METHOD_LOCK_SHA256.txt`.
- The generic preceding numerical engine is vendored under
  `reference_implementation/`; no source is loaded from a sibling directory.
- Every data, calibration-split, CV-fold, and independent-test seed is an
  explicit deterministic function of the master seed, scenario, and
  replication. The committed production ledger contains all 7,000 rows.
- Production, grid, dedicated-smoke, development, and retired pre-lock seed
  families are enumerated and checked for exact non-overlap.
- A calibration failure retains the prespecified `kappa=1` default and a flag;
  no failed value is silently replaced. The production run used no defaults.
- Every denominator in the paired relative test-deviance estimand must be
  finite and strictly positive. Analysis stops on a violation: no replication
  is discarded and no denominator floor or replacement is applied. The
  realised denominator ranges and relative-contrast second moments are retained
  in `results/relative_deviance_denominator_audit.csv`.
- The nonlinear logistic and Poisson active-count rows are labelled BIC-like
  for compactness but remain exploratory analogues outside the Gaussian
  BIC-like proposition.

See `METHODS_AND_OUTPUTS.md` for the formal estimands and output schema and
`PROVENANCE.md` for the development/confirmation separation and packaging
record.

## Source files

| File | Purpose |
| --- | --- |
| `Simulation_Study.R` | Main validation, simulation, analysis, and verification pipeline |
| `Simulation_Config.R` | Locked scenarios, grids, seeds, margins, and solver settings |
| `PanIC_CF_Functions.R` | Revised weight, same-path original sensitivity, CV-min, and CV-1SE |
| `Run_Confirmatory_Simulations.R` | Parallel seven-setting runner and exact seed ledger |
| `Simulation_Results_Analysis.R` | Monte Carlo summaries, joint decision, denominator audit, and generated production tables |
| `Results_Verification.R` | Independent checks of full replication-level main results and denominator audit |
| `Grid_Sensitivity_Study.R` | Fresh-seed common-random-number grid experiment |
| `Grid_Sensitivity_Pipeline.R` | Grid validation, execution, and verification entry point |
| `Grid_Results_Verification.R` | Full replication-level grid checks |
| `Implementation_Validation.R` | Deterministic method, equivalence, and seed checks |
| `Method_Lock_Verification.R` | SHA-256 verification of every locked source |
| `Repository_Verification.R` | Verification of the compact committed release |
| `reference_implementation/` | Frozen generic numerical engine inherited by PanIC-CF |
| `Boundary_Subsampling_Study.R` | Scalar boundary and repeated-subsampling study |
| `Solver_Validation.R` | Direct-fit and fixed-shape Gamma diagnostics |
| `Runtime_Benchmark.R` | One-worker computational benchmark |
| `Render_Runtime_Figure.R` | Runtime figure renderer using archived timings |
| `Manuscript_Table_Rendering.R` | Pure deterministic renderers for all seven retained table mirrors |
| `Render_Manuscript_Tables.R` | Command-line table-rendering entry point |
| `Manuscript_Table_Tools.R` | Exact refresh and verification of marker-delimited tables embedded in `manuscript/main.tex` |
| `manuscript/` | Self-contained manuscript, reviewer response, bibliography, and figures |

## Reference outputs

The production run contains 7,000 data sets, 35,000 primary method rows, and
70,000 calibration-direction rows. It recorded no method, path, calibration,
or CV failure; no default use; no solver warning; no multiplier endpoint; and
no PanIC-CF or CV selected-radius endpoint. One of 70,000 calibration targets
was projected to the allowed radius interval. All 21 deterministic checks and
the independent main and grid result verifiers passed.

| Manuscript item | Retained files |
| --- | --- |
| Confirmatory decision | `results/confirmatory_decision.csv`, `results/scenario_primary_estimands.csv`, `results/relative_deviance_denominator_audit.csv`, `results/table_confirmatory_decision.tex` |
| Primary support and performance | `results/simulation_summary.csv`, `results/paired_method_contrasts.csv`, `results/table_primary_support.tex`, `results/table_primary_performance.tex` |
| Calibration diagnostics | `results/calibration_target_summary.csv`, `results/diagnostic_summary.csv`, `results/table_calibration.tex` |
| Grid sensitivity | `results/grid_sensitivity_summary.csv`, `results/grid_sensitivity_paired_contrasts.csv`, `results/grid_sensitivity_method_contrasts.csv`, `results/table_grid_sensitivity.tex` |
| Boundary and subsampling | `results/boundary_summary.csv`, `results/subsampling_summary.csv`, `results/table_boundary.tex`, `results/table_subsampling.tex`, `results/figure_boundary_subsampling.pdf` |
| Runtime | `results/runtime_summary.csv`, `results/runtime_environment.txt`, `results/figure_runtime.pdf` |
| Solver and Gamma checks | `results/solver_validation.csv`, `results/gamma_validation.csv` |

The `results/table_*.tex` files in this table are reproducibility mirrors of
the seven blocks embedded in `manuscript/main.tex`; the manuscript does not
load them at compile time.

## Citation and contact

Citation metadata are provided in `CITATION.cff`; no DOI has been assigned in
this repository.

Hien Duy Nguyen: h.nguyen5@latrobe.edu.au
