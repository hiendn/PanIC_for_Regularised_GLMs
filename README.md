# Consistent information criteria for regularised generalised linear models

R code and simulation results for the paper by Qingyuan Zhang and Hien Duy Nguyen.

## Requirements

The reference results use R 4.4.0, `glmnet` 4.1-10 and `Matrix` 1.7-0 on
`aarch64-apple-darwin20`. Install `glmnet` if needed:

```r
install.packages("glmnet")
```

## Reproduce and verify

Run from the repository root in a fresh checkout:

```sh
Rscript Implementation_Validation.R
Rscript Repository_Verification.R
```

The first command checks the simulation design, fitting rules and random-number
streams. The second checks the saved summaries, tables, figures and manuscript
against the retained results. Both commands stop if a check fails.

To rerun the primary simulation:

```sh
Rscript Simulation_Study.R --n-rep=1000 --cores=4 \
  --output-dir=production_results --overwrite
```

The primary study compares PanIC-CF, ordinary five-fold cross-validation and a
BIC-like comparator. It uses five data-generating designs at sample sizes 500
and 1,000, with 1,000 replications per setting. The master seed is
`2146092101`.

The additional grid, boundary, subsampling and runtime studies can be rerun
with:

```sh
Rscript Grid_Sensitivity_Pipeline.R --n-rep=500 --cores=4 \
  --output-dir=production_results
Rscript Boundary_Subsampling_Study.R --output-dir=production_results
Rscript Runtime_Benchmark.R --output-dir=production_results
```

Full runs write their detailed outputs to `production_results/`. The smaller
`results/` directory contains the summaries and checks used by the manuscript.

## Source files

| File | Purpose |
| --- | --- |
| `Simulation_Study.R` | Main entry point for the primary simulation |
| `Simulation_Config.R` | Designs, seeds, grids and fitting settings |
| `PanIC_CF_Functions.R` | Data generation, model fitting and selection methods |
| `Run_Confirmatory_Simulations.R` | Simulation runner |
| `Simulation_Results_Analysis.R` | Summary tables and paired comparisons |
| `Results_Verification.R` | Checks for the primary simulation results |
| `Grid_Sensitivity_Pipeline.R` | Grid-sensitivity study |
| `Boundary_Subsampling_Study.R` | Boundary and subsampling illustrations |
| `Runtime_Benchmark.R` | Runtime comparison |
| `Implementation_Validation.R` | Implementation checks |
| `Repository_Verification.R` | Verification of the retained repository outputs |
| `manuscript/` | Manuscript and response to reviewers |

## Reference outputs

The numerical implementation used for the reference outputs is identified by
[source commit d7bd8d3](https://github.com/hiendn/PanIC_for_Regularised_GLMs/commit/d7bd8d3ece811a28eca4fe52487d05250e7b7f5b).
The reference run contains 10,000 simulated data sets and 30,000 primary
method rows. All checks in `results/repository_verification.csv` passed.

| Manuscript result | File |
| --- | --- |
| Support recovery | `results/table_primary_support.tex` |
| Selected radii and test performance | `results/table_primary_performance.tex` |
| PanIC-CF calibration | `results/table_calibration.tex` |
| Radius-grid sensitivity | `results/table_grid_sensitivity.tex` |
| Boundary illustration | `results/table_boundary.tex` |
| Subsampling illustration | `results/table_subsampling.tex` |

The runtime and boundary/subsampling figures are also retained under
`results/`. Citation metadata are provided in `CITATION.cff`.

## Contact

Hien Duy Nguyen: h.nguyen5@latrobe.edu.au
