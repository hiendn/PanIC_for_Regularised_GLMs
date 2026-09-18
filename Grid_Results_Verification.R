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
} else {
  normalizePath("Grid_Results_Verification.R")
}
candidate_dir <- dirname(script_path)
source(file.path(candidate_dir, "Method_Lock_Verification.R"))
verify_method_lock(candidate_dir)
source(file.path(candidate_dir, "Simulation_Config.R"))
output_name <- arg_value("--output-dir", CONFIG$default_output_dir)
results_dir <- if (grepl("^/", output_name)) {
  output_name
} else {
  file.path(candidate_dir, output_name)
}

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
read_result <- function(name) {
  path <- file.path(results_dir, name)
  assert(file.exists(path), paste("Missing grid result:", name))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

grid_configuration <- read_result("grid_sensitivity_configuration.csv")
n_rep <- as.integer(grid_configuration$value[
  grid_configuration$item == "replications"
])
seed_role <- grid_configuration$value[
  grid_configuration$item == "seed_role"
]
assert(length(seed_role) == 1L &&
         seed_role %in% c("production", "dedicated-smoke"),
       "Grid seed role is missing or invalid")
expected_grid_seed <- if (seed_role == "production") {
  CONFIG$grid_sensitivity_master_seed
} else {
  CONFIG$smoke_grid_sensitivity_master_seed
}
expected_main_seed <- if (seed_role == "production") {
  CONFIG$master_seed
} else {
  CONFIG$smoke_master_seed
}
assert(as.integer(grid_configuration$value[
  grid_configuration$item == "sensitivity_master_seed"
]) == expected_grid_seed,
"Grid-sensitivity master seed is not locked")
assert(as.integer(grid_configuration$value[
  grid_configuration$item == "confirmatory_master_seed"
]) == expected_main_seed, "Main master seed changed in the grid archive")
assert(as.integer(grid_configuration$value[
  grid_configuration$item == "production_sensitivity_master_seed"
]) == CONFIG$grid_sensitivity_master_seed,
"Recorded production grid seed changed")
assert(as.integer(grid_configuration$value[
  grid_configuration$item == "production_confirmatory_master_seed"
]) == CONFIG$master_seed, "Recorded production main seed changed")
if (seed_role == "production") {
  assert(n_rep == CONFIG$grid_sensitivity_replications,
         "Production grid archive has the wrong replication count")
}

method_order <- c(
  "PanIC-CF", CONFIG$original_sensitivity_method, "BIC-like",
  CONFIG$primary_cv_method, CONFIG$secondary_cv_method
)
calibration_by_grid <- list()
for (m in CONFIG$grid_sensitivity_points) {
  primary <- read_result(paste0("grid_sensitivity_m", m, "_primary.csv"))
  diagnostic <- read_result(
    paste0("grid_sensitivity_m", m, "_diagnostics.csv")
  )
  calibration <- read_result(
    paste0("grid_sensitivity_m", m, "_calibration_rows.csv")
  )
  calibration_by_grid[[as.character(m)]] <- calibration
  assert(nrow(primary) == 5L * n_rep,
         paste0("m=", m, ": primary row count mismatch"))
  assert(nrow(diagnostic) == n_rep,
         paste0("m=", m, ": diagnostic row count mismatch"))
  assert(nrow(calibration) == 10L * n_rep,
         paste0("m=", m, ": calibration row count mismatch"))
  assert(all(vapply(
    split(primary$method, primary$replication),
    function(methods) setequal(methods, method_order), logical(1)
  )), paste0("m=", m, ": method set mismatch"))
  cv_min <- primary[primary$method == CONFIG$primary_cv_method, ]
  cv_one_se <- primary[primary$method == CONFIG$secondary_cv_method, ]
  ok <- cv_min$method_failed == 0L & cv_one_se$method_failed == 0L
  assert(all(
    cv_one_se$selected_grid_radius[ok] <=
      cv_min$selected_grid_radius[ok] + 1e-14
  ), paste0("m=", m, ": CV-1SE radius exceeds CV-min"))
  expected_weight <- sqrt(log(log(
    calibration$n_validation + exp(exp(1))
  )))
  expected_original <- log(log(
    calibration$n_validation + exp(exp(1))
  ))
  target_ok <- is.finite(calibration$raw_cross_signed_radius)
  assert(max(abs(
    calibration$weight[target_ok] - expected_weight[target_ok]
  )) < 1e-12, paste0("m=", m, ": revised weight mismatch"))
  assert(max(abs(
    calibration$original_weight[target_ok] - expected_original[target_ok]
  )) < 1e-12, paste0("m=", m, ": original weight mismatch"))
}

reference <- calibration_by_grid[["121"]]
keys <- c("replication", "split_repeat", "direction")
for (m in c("61", "241")) {
  candidate <- calibration_by_grid[[m]]
  pair <- merge(reference, candidate, by = keys,
                suffixes = c("_reference", "_candidate"))
  assert(nrow(pair) == 10L * n_rep,
         paste0("m=", m, ": common-random-number pairing incomplete"))
  for (field in c(
    "training_seed", "split_seed", "test_seed", "raw_cross_signed_radius",
    "weight", "original_weight"
  )) {
    lhs <- pair[[paste0(field, "_reference")]]
    rhs <- pair[[paste0(field, "_candidate")]]
    assert(isTRUE(all.equal(lhs, rhs, tolerance = 1e-12)),
           paste0("m=", m, ": common quantity changed: ", field))
  }
}

summary <- read_result("grid_sensitivity_summary.csv")
grid_contrasts <- read_result("grid_sensitivity_paired_contrasts.csv")
method_contrasts <- read_result("grid_sensitivity_method_contrasts.csv")
assert(nrow(summary) == 15L,
       "Grid summary must contain five methods on three grids")
assert(nrow(grid_contrasts) == 15L,
       "Grid paired table must contain five methods and three grid pairs")
assert(nrow(method_contrasts) == 12L,
       "Grid method table must contain four comparisons on three grids")
assert(all(summary$attempted_replications == n_rep),
       "Grid summary replication count mismatch")
assert(file.exists(file.path(
  results_dir, "manuscript_generated", "table_grid_sensitivity.tex"
)), "Missing grid-sensitivity LaTeX table")
assert(file.exists(file.path(results_dir, "grid_locked_configuration.rds")),
       "Missing locked grid configuration")
grid_table_text <- paste(readLines(
  file.path(results_dir, "manuscript_generated", "table_grid_sensitivity.tex"),
  warn = FALSE
), collapse = "\n")
assert(!grepl("BIC-active \\(exploratory\\)", grid_table_text) &&
         grepl("BIC-like", grid_table_text),
       "Reader-facing grid BIC-like label was not normalized")
assert(!grepl(
  "second-confirmation|revised PanIC-CF", grid_table_text,
  ignore.case = TRUE
), "The grid table exposes internal process wording")

cat("PASS: confirmatory grid-sensitivity results are verified.\n")
