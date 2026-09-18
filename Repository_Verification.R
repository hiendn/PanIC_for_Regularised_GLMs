#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[[1L]])
}
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]))
} else normalizePath("Repository_Verification.R")
repository_dir <- dirname(script_path)
output_name <- arg_value("--output-dir", "results")
results_dir <- if (grepl("^/", output_name)) output_name else
  file.path(repository_dir, output_name)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
read_result <- function(name) {
  path <- file.path(results_dir, name)
  assert(file.exists(path), paste("Missing retained result:", name))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

scenarios <- read_result("scenario_registry.csv")
configuration <- read_result("configuration.csv")
summary <- read_result("simulation_summary.csv")
paired <- read_result("paired_panic_cf_vs_cv.csv")
diagnostics <- read_result("diagnostic_summary.csv")
targets <- read_result("calibration_target_summary.csv")
grid_configuration <- read_result("grid_sensitivity_configuration.csv")
grid_summary <- read_result("grid_sensitivity_summary.csv")
grid_paired <- read_result("grid_sensitivity_paired_contrasts.csv")
runtime <- read_result("runtime_summary.csv")
checks <- read_result("verification_checks.csv")

assert(nrow(scenarios) == 7L, "Scenario registry must contain seven settings")
assert(nrow(summary) == 21L, "Primary summary must contain three methods in seven settings")
assert(nrow(paired) == 7L, "Paired PanIC-CF/CV table must contain seven settings")
assert(nrow(diagnostics) == 7L && nrow(targets) == 7L,
       "Diagnostic and calibration summaries must contain seven settings")
assert(all(summary$attempted_replications == 500L) &&
         all(summary$failed_replications == 0L),
       "Primary summaries do not record 500 successful attempts per method")
assert(all(paired$paired_replications == 500L),
       "A paired comparison does not contain all 500 replications")

expected_support <- c(-0.472, -0.428, 0.190, 0.090, -0.426, 0.050, -0.360)
paired <- paired[match(scenarios$scenario_id, paired$scenario_id), ]
assert(max(abs(paired$wrong_difference - expected_support)) < 1e-12,
       "A retained total-support-error contrast has changed")
## The unrounded maximum is 1.1321908e-4, reported as 1.13e-4.
assert(all(abs(paired$test_deviance_difference) <= 1.14e-4) &&
         all(abs(paired$test_deviance_difference) <
               paired$test_deviance_difference_mcse),
       "The retained predictive-comparison qualification is not satisfied")

zero_diagnostic_columns <- c(
  "full_path_failures", "calibration_failures",
  "calibration_default_uses", "cv_failures",
  "raw_targets_below_grid", "raw_targets_above_grid", "projections",
  "kappa_lower_endpoint_selections", "kappa_upper_endpoint_selections",
  "solver_warnings"
)
assert(all(vapply(diagnostics[zero_diagnostic_columns], function(x) all(x == 0),
                  logical(1))),
       "A retained confirmatory diagnostic is nonzero")
assert(all(targets$complete_replications == 500L) &&
         max(abs(targets$raw_target_mean -
                 targets$projected_target_mean)) < 1e-12,
       "Calibration-target summaries are incomplete or unexpectedly projected")

assert(nrow(grid_summary) == 9L && nrow(grid_paired) == 9L,
       "Grid-sensitivity summaries must each contain nine rows")
assert(all(grid_summary$attempted_replications == 500L) &&
         all(grid_summary$method_failures == 0L),
       "Grid-sensitivity reference results are incomplete")
assert(as.integer(grid_configuration$value[
  grid_configuration$item == "sensitivity_master_seed"
]) == 2036091802L, "Grid-sensitivity seed has changed")
assert(nrow(runtime) == 9L && all(runtime$attempted_replications == 20L) &&
         all(runtime$failures == 0L),
       "Runtime summary must contain nine complete 20-run cells")
assert(nrow(checks) == 14L && all(checks$passed == 1L),
       "One or more retained deterministic checks failed")

boundary <- read_result("boundary_summary.csv")
subsampling <- read_result("subsampling_summary.csv")
solver <- read_result("solver_validation.csv")
gamma <- read_result("gamma_validation.csv")
assert(nrow(boundary) == 2L && all(boundary$replications == 40000L),
       "Boundary summary is incomplete")
assert(nrow(subsampling) == 3L && all(subsampling$random_subsamples == 10000L),
       "Subsampling summary is incomplete")
assert(nrow(solver) == 9L, "Direct solver-validation summary must contain nine rows")
assert(nrow(gamma) == 1L && gamma$warning_count == 1L &&
         gamma$glmnet_jerr == 0L,
       "The fixed-shape Gamma diagnostic no longer matches its disclosed scope")

manuscript_files <- c(
  "table_primary_support.tex", "table_primary_performance.tex",
  "table_calibration.tex", "table_grid_sensitivity.tex",
  "table_boundary.tex", "table_subsampling.tex",
  "figure_runtime.pdf", "figure_boundary_subsampling.pdf"
)
paths <- file.path(results_dir, manuscript_files)
assert(all(file.exists(paths)) && all(file.info(paths)$size > 0),
       "A retained manuscript-facing table or figure is missing")

cat("PASS: compact repository outputs and numerical claims are verified.\n")
