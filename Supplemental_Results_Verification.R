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
} else normalizePath("Supplemental_Results_Verification.R")
production_dir <- dirname(script_path)
source(file.path(production_dir, "Simulation_Config.R"))
output_name <- arg_value("--output-dir", "results")
results_dir <- if (grepl("^/", output_name)) output_name else
  file.path(production_dir, output_name)
assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

runtime <- read.csv(file.path(results_dir, "runtime_raw.csv"),
                    stringsAsFactors = FALSE)
runtime_summary <- read.csv(file.path(results_dir, "runtime_summary.csv"),
                            stringsAsFactors = FALSE)
assert(nrow(runtime) == 3L * 3L * CONFIG$runtime_replications,
       "Runtime archive has an unexpected row count")
runtime_cell_counts <- table(interaction(
  runtime$family, runtime$method, drop = TRUE
))
assert(length(runtime_cell_counts) == 9L &&
         all(runtime_cell_counts == CONFIG$runtime_replications),
       "A runtime family/method cell is incomplete")
assert(all(runtime$failed == 0L), "A runtime replication failed")
expected_paths <- c(
  `Revised PanIC-CF` = 11L, `5-fold CV` = 6L, `BIC-like` = 1L,
  `Exploratory active-count` = 1L
)
for (method in names(expected_paths)) {
  dat <- runtime[runtime$method == method, ]
  assert(all(dat$paths_fitted == expected_paths[[method]]),
         paste0(method, " has an incorrect path count"))
}
panic_runtime <- runtime[runtime$method == "Revised PanIC-CF", ]
assert(all(panic_runtime$unregularised_sign_fits == 10L) &&
         all(panic_runtime$unregularised_validation_fits == 10L),
       "Revised PanIC-CF runtime does not include all pilot GLMs")
recomputed <- aggregate(
  elapsed_seconds ~ family + n + method, runtime,
  function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)),
                median = median(x), minimum = min(x), maximum = max(x))
)
stats <- recomputed$elapsed_seconds
if (is.list(stats)) stats <- do.call(rbind, stats)
stats <- as.data.frame(stats)
recomputed <- cbind(recomputed[c("family", "n", "method")], stats)
check <- merge(
  recomputed, runtime_summary,
  by = c("family", "n", "method"), suffixes = c("_check", "_saved")
)
for (metric in c("mean", "se", "median", "minimum", "maximum")) {
  assert(max(abs(check[[paste0(metric, "_check")]] -
                 check[[paste0(metric, "_saved")]])) < 1e-12,
         paste0("Runtime summary mismatch for ", metric))
}
assert(file.exists(file.path(
  results_dir, "figure_runtime.pdf"
)), "Revised runtime figure is missing")
assert(file.exists(file.path(results_dir, "runtime_environment.txt")),
       "Runtime environment record is missing")

grid_sizes <- c(61L, 121L, 241L)
calibration_sets <- list()
for (m in grid_sizes) {
  primary <- read.csv(file.path(
    results_dir, paste0("grid_sensitivity_m", m, "_primary.csv")
  ), stringsAsFactors = FALSE)
  diagnostic <- read.csv(file.path(
    results_dir, paste0("grid_sensitivity_m", m, "_diagnostics.csv")
  ), stringsAsFactors = FALSE)
  calibration <- read.csv(file.path(
    results_dir, paste0("grid_sensitivity_m", m,
                       "_calibration_rows.csv")
  ), stringsAsFactors = FALSE)
  calibration_sets[[as.character(m)]] <- calibration
  assert(nrow(primary) == 3L * 500L,
         paste0("Grid m=", m, " primary row count is incomplete"))
  assert(nrow(diagnostic) == 500L,
         paste0("Grid m=", m, " diagnostics are incomplete"))
  assert(nrow(calibration) == 10L * 500L,
         paste0("Grid m=", m, " calibration rows are incomplete"))
  assert(all(diagnostic$full_path_failed == 0L),
         paste0("Grid m=", m, " has a full-path failure"))
  assert(all(diagnostic$calibration_failed == 0L),
         paste0("Grid m=", m, " has a calibration failure"))
  assert(all(diagnostic$default_used == 0L),
         paste0("Grid m=", m, " used the calibration default"))
  assert(all(diagnostic$cv_failed == 0L),
         paste0("Grid m=", m, " has a CV failure"))
  assert(max(diagnostic$maximum_active_radius_error) <=
           CONFIG$interpolation_radius_tolerance,
         paste0("Grid m=", m, " exceeds the interpolation tolerance"))
}

key <- c("replication", "split_repeat", "direction")
reference <- calibration_sets[["121"]]
for (m in c("61", "241")) {
  comparison <- merge(
    reference[c(key, "training_seed", "split_seed", "test_seed",
                "raw_cross_signed_radius")],
    calibration_sets[[m]][c(
      key, "training_seed", "split_seed", "test_seed",
      "raw_cross_signed_radius"
    )],
    by = key, suffixes = c("_reference", "_candidate")
  )
  assert(nrow(comparison) == 5000L,
         paste0("Grid m=", m, " common-random-number join is incomplete"))
  assert(all(comparison$training_seed_reference ==
             comparison$training_seed_candidate) &&
         all(comparison$split_seed_reference ==
             comparison$split_seed_candidate) &&
         all(comparison$test_seed_reference ==
             comparison$test_seed_candidate),
         paste0("Grid m=", m, " does not reuse the same random streams"))
  assert(max(abs(comparison$raw_cross_signed_radius_reference -
                 comparison$raw_cross_signed_radius_candidate)) < 1e-12,
         paste0("Grid m=", m, " raw targets differ under common data"))
}

grid_configuration <- read.csv(file.path(
  results_dir, "grid_sensitivity_configuration.csv"
), stringsAsFactors = FALSE)
sensitivity_seed <- as.integer(grid_configuration$value[
  grid_configuration$item == "sensitivity_master_seed"
])
confirmatory_seed <- as.integer(grid_configuration$value[
  grid_configuration$item == "confirmatory_master_seed"
])
assert(sensitivity_seed == 2036091802L &&
         confirmatory_seed == CONFIG$master_seed &&
         sensitivity_seed != confirmatory_seed,
       "Sensitivity and confirmatory seed streams are not as prespecified")
grid_summary <- read.csv(file.path(
  results_dir, "grid_sensitivity_summary.csv"
), stringsAsFactors = FALSE)
paired <- read.csv(file.path(
  results_dir, "grid_sensitivity_paired_contrasts.csv"
), stringsAsFactors = FALSE)
assert(nrow(grid_summary) == 9L, "Grid summary must have nine rows")
assert(nrow(paired) == 9L, "Grid paired-contrast table must have nine rows")
assert(all(paired$paired_replications == 500L),
       "A grid paired contrast does not use all replications")
assert(file.exists(file.path(
  results_dir, "table_grid_sensitivity.tex"
)), "Revised grid-sensitivity table is missing")

cat("All runtime and grid-sensitivity integrity checks passed.\n")
