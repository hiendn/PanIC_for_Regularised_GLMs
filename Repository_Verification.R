#!/usr/bin/env Rscript

## Verify the compact, committed PanIC-CF release without requiring the bulky
## replication-level production files.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]))
} else {
  normalizePath("Repository_Verification.R")
}
repository_dir <- dirname(script_path)
results_dir <- file.path(repository_dir, "results")

source(file.path(repository_dir, "Method_Lock_Verification.R"))
verify_method_lock(repository_dir)
source(file.path(repository_dir, "Simulation_Config.R"))
source(file.path(repository_dir, "Manuscript_Table_Rendering.R"))
source(file.path(repository_dir, "Manuscript_Table_Tools.R"))

checks <- list()
record_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = as.integer(isTRUE(passed)),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}
assert_file <- function(name) {
  path <- file.path(results_dir, name)
  if (!file.exists(path)) stop("Missing committed result: ", name, call. = FALSE)
  path
}
read_result <- function(name) {
  read.csv(assert_file(name), stringsAsFactors = FALSE, check.names = FALSE)
}
sha256_file <- function(path) {
  if (nzchar(Sys.which("shasum"))) {
    output <- system2(
      "shasum", c("-a", "256", shQuote(path)), stdout = TRUE
    )
  } else if (nzchar(Sys.which("sha256sum"))) {
    output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  } else {
    stop("Neither shasum nor sha256sum is available", call. = FALSE)
  }
  substr(output[[1L]], 1L, 64L)
}
close_enough <- function(x, y, tolerance = 5e-13) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance,
                   check.attributes = FALSE))
}
assessed_methods <- c("PanIC-CF", "BIC-like", "CV")
obsolete_method_pattern <- paste(
  c(
    "PanIC-CF-original", "PanIC-CF-loglog", "CV-min", "CV-1SE",
    "BIC-active"
  ),
  collapse = "|"
)
contains_obsolete_method_label <- function(data) {
  character_columns <- vapply(data, is.character, logical(1))
  if (!any(character_columns)) return(FALSE)
  values <- unlist(data[character_columns], use.names = FALSE)
  any(grepl(obsolete_method_pattern, values))
}

## Check the compact release checksum manifest first.
checksum_path <- assert_file("RESULTS_SHA256.txt")
checksum_lines <- readLines(checksum_path, warn = FALSE)
checksum_entries <- grep("^[0-9a-f]{64}  ", checksum_lines, value = TRUE)
checksum_expected <- substr(checksum_entries, 1L, 64L)
checksum_relative <- substring(checksum_entries, 67L)
expected_checksum_basenames <- c(
  "production_file_manifest.csv",
  "implementation_validation.csv",
  "configuration.csv",
  "scenario_registry.csv",
  "locked_configuration.rds",
  "environment.txt",
  "sessionInfo.txt",
  "seed_ledger.csv",
  "simulation_summary.csv",
  "paired_method_contrasts.csv",
  "paired_panic_cf_vs_cv.csv",
  "scenario_primary_estimands.csv",
  "confirmatory_decision.csv",
  "prediction_margin_sensitivity.csv",
  "relative_deviance_denominator_audit.csv",
  "calibration_target_summary.csv",
  "diagnostic_summary.csv",
  "grid_sensitivity_configuration.csv",
  "grid_locked_configuration.rds",
  "grid_sensitivity_summary.csv",
  "grid_sensitivity_paired_contrasts.csv",
  "grid_sensitivity_method_contrasts.csv",
  "runtime_summary.csv",
  "runtime_environment.txt",
  "boundary_summary.csv",
  "subsampling_summary.csv",
  "solver_validation.csv",
  "gamma_validation.csv",
  PANIC_RENDERED_TABLE_FILES,
  "figure_runtime.pdf",
  "figure_boundary_subsampling.pdf"
)
expected_checksum_relative <- file.path(
  "results", expected_checksum_basenames
)
checksum_inventory_ok <-
  !anyDuplicated(checksum_relative) &&
  identical(sort(checksum_relative), sort(expected_checksum_relative))
checksum_files <- file.path(repository_dir, checksum_relative)
checksum_present <- checksum_inventory_ok && all(file.exists(checksum_files))
checksum_actual <- if (checksum_present) {
  vapply(checksum_files, sha256_file, character(1))
} else {
  rep(NA_character_, length(checksum_expected))
}
record_check(
  "compact result checksum inventory and digests",
  checksum_present && identical(unname(checksum_actual), checksum_expected),
  paste0(
    "entries=", length(checksum_entries),
    ";missing=", paste(
      setdiff(expected_checksum_relative, checksum_relative), collapse = ","
    ),
    ";extra=", paste(
      setdiff(checksum_relative, expected_checksum_relative), collapse = ","
    )
  )
)
result_file_names <- list.files(results_dir)
obsolete_result_files <- unique(c(grep(
  "cv[_-]?min|cv[_-]?1se|loglog|panic[^/]*original|bic[_-]?active",
  result_file_names, value = TRUE, ignore.case = TRUE
), intersect(result_file_names, "table_confirmatory_decision.tex")))
record_check(
  "obsolete compact-result filenames are absent",
  !length(obsolete_result_files),
  if (length(obsolete_result_files)) {
    paste(obsolete_result_files, collapse = ",")
  } else {
    "none"
  }
)

implementation <- read_result("implementation_validation.csv")
record_check(
  "all deterministic implementation checks passed",
  nrow(implementation) == 19L && all(implementation$passed == 1L),
  paste0("passed=", sum(implementation$passed), "/", nrow(implementation))
)

configuration <- read_result("configuration.csv")
config_value <- function(item) {
  value <- configuration$value[configuration$item == item]
  if (length(value) != 1L) stop("Configuration item is not unique: ", item)
  value
}
record_check(
  "production configuration and diagnostic thresholds",
  as.integer(config_value("master_seed")) == 2146092101L &&
    as.integer(config_value("production_master_seed")) == 2146092101L &&
    as.integer(config_value("replications_requested")) == 1000L &&
    as.integer(config_value("independent_test_size")) == 2000L &&
    as.numeric(config_value("prediction_noninferiority_margin")) == 0.001 &&
    config_value("cv_primary_rule") ==
      "minimum mean five-fold validation loss on radius grid" &&
    config_value("support_definition") ==
      "beta_hat != 0 (literal fitted nonzero)" &&
    !contains_obsolete_method_label(configuration),
  paste0("seed=", config_value("master_seed"),
         ";replications=", config_value("replications_requested"))
)

registry <- read_result("scenario_registry.csv")
record_check(
  "ten simulation scenarios are present",
  nrow(registry) == nrow(SCENARIOS) &&
    identical(registry$scenario_id, SCENARIOS$scenario_id) &&
    identical(registry$family, SCENARIOS$family) &&
    identical(as.integer(registry$n), as.integer(SCENARIOS$n)) &&
    close_enough(registry$rho, SCENARIOS$rho),
  paste(registry$scenario_id, collapse = ",")
)

seed_ledger <- read_result("seed_ledger.csv")
seed_columns <- c(
  "scenario_id", "replication", "training_seed", "cv_fold_seed", "test_seed",
  paste0("calibration_split_seed_", seq_len(CONFIG$calibration_half_splits))
)
seed_schema_ok <- all(seed_columns %in% names(seed_ledger))
seed_values_ok <- seed_schema_ok &&
  nrow(seed_ledger) == nrow(SCENARIOS) * CONFIG$n_rep
if (seed_values_ok) {
  for (i in seq_len(nrow(seed_ledger))) {
    scenario_index <- match(seed_ledger$scenario_id[i], SCENARIOS$scenario_id)
    expected <- seed_streams(scenario_index, seed_ledger$replication[i])
    observed_splits <- as.integer(unlist(seed_ledger[i, paste0(
      "calibration_split_seed_", seq_len(CONFIG$calibration_half_splits)
    )], use.names = FALSE))
    if (seed_ledger$training_seed[i] != expected$training ||
        !identical(observed_splits, expected$calibration_splits) ||
        seed_ledger$cv_fold_seed[i] != expected$cv_folds ||
        seed_ledger$test_seed[i] != expected$test) {
      seed_values_ok <- FALSE
      break
    }
  }
}
record_check(
  "complete deterministic production seed ledger",
  seed_values_ok && !anyDuplicated(seed_ledger[c("scenario_id", "replication")]),
  paste0("rows=", nrow(seed_ledger))
)

summary <- read_result("simulation_summary.csv")
method_sets_ok <- all(vapply(split(summary, summary$scenario_id), function(dat) {
  setequal(dat$method, assessed_methods)
}, logical(1)))
record_check(
  "primary summary is complete",
  nrow(summary) == length(assessed_methods) * nrow(SCENARIOS) &&
    method_sets_ok &&
    all(summary$attempted_replications == CONFIG$n_rep) &&
    all(summary$failed_replications == 0L) &&
    all(summary$method %in% assessed_methods) &&
    !contains_obsolete_method_label(summary),
  paste0("rows=", nrow(summary), ";failures=",
         sum(summary$failed_replications))
)

paired <- read_result("paired_method_contrasts.csv")
primary_pair <- read_result("paired_panic_cf_vs_cv.csv")
expected_pair_keys <- c(
  "PanIC-CF\rCV", "PanIC-CF\rBIC-like", "BIC-like\rCV"
)
paired_sets_ok <- all(vapply(split(paired, paired$scenario_id), function(dat) {
  keys <- paste(dat$lhs_method, dat$rhs_method, sep = "\r")
  setequal(keys, expected_pair_keys)
}, logical(1)))
record_check(
  "compact paired summaries use the assessed methods",
  nrow(paired) == length(expected_pair_keys) * nrow(SCENARIOS) &&
    paired_sets_ok &&
    all(paired$paired_replications == CONFIG$n_rep) &&
    nrow(primary_pair) == nrow(SCENARIOS) &&
    identical(primary_pair$scenario_id, SCENARIOS$scenario_id) &&
    all(primary_pair$lhs_method == "PanIC-CF") &&
    all(primary_pair$rhs_method == "CV") &&
    all(primary_pair$paired_replications == CONFIG$n_rep) &&
    !contains_obsolete_method_label(paired) &&
    !contains_obsolete_method_label(primary_pair),
  paste0("all_pairs=", nrow(paired), ";primary_pairs=", nrow(primary_pair))
)

scenario <- read_result("scenario_primary_estimands.csv")
decision <- read_result("confirmatory_decision.csv")
z <- qnorm(0.95)
support_estimate <- mean(scenario$wrong_difference)
support_mcse <- sqrt(sum(scenario$wrong_difference_mcse^2)) / nrow(scenario)
support_upper <- support_estimate + z * support_mcse
prediction_estimate <- mean(scenario$relative_test_deviance_difference)
prediction_mcse <- sqrt(sum(
  scenario$relative_test_deviance_difference_mcse^2
)) / nrow(scenario)
prediction_upper <- prediction_estimate + z * prediction_mcse
expected_complete_pairing <- all(
  scenario$paired_replications == CONFIG$n_rep
)
expected_support_pass <- is.finite(support_upper) && support_upper < 0
expected_prediction_pass <- is.finite(prediction_upper) &&
  prediction_upper < CONFIG$prediction_noninferiority_margin
expected_joint_flag <- expected_complete_pairing && expected_support_pass &&
  expected_prediction_pass
decision_values_ok <-
  nrow(scenario) == nrow(SCENARIOS) && nrow(decision) == 1L &&
  all(scenario$paired_replications == CONFIG$n_rep) &&
  decision$scenarios == nrow(SCENARIOS) &&
  decision$replications_per_scenario == CONFIG$n_rep &&
  close_enough(decision$support_estimate, support_estimate) &&
  close_enough(decision$support_mcse, support_mcse) &&
  close_enough(decision$support_upper_one_sided_95, support_upper) &&
  close_enough(decision$prediction_estimate, prediction_estimate) &&
  close_enough(decision$prediction_mcse, prediction_mcse) &&
  close_enough(decision$prediction_upper_one_sided_95, prediction_upper) &&
  decision$complete_pairing == as.integer(expected_complete_pairing) &&
  decision$support_superiority_pass == as.integer(expected_support_pass) &&
  decision$prediction_noninferiority_pass ==
    as.integer(expected_prediction_pass) &&
  decision$joint_claim_pass == as.integer(expected_joint_flag)
record_check(
  "endpoint summaries and flags reconstruct from scenario rows",
  decision_values_ok,
  sprintf("support_upper=%.9f;prediction_upper=%.9f",
          support_upper, prediction_upper)
)

denominator_audit <- read_result("relative_deviance_denominator_audit.csv")
denominator_audit_columns <- c(
  "scenario_id", "family", "n", "rho", "paired_replications",
  "nonfinite_cv_test_deviances",
  "nonpositive_cv_test_deviances",
  "minimum_cv_test_deviance",
  "maximum_cv_test_deviance",
  "nonfinite_relative_contrasts",
  "relative_contrast_mean",
  "relative_contrast_variance",
  "relative_contrast_second_moment",
  "maximum_absolute_relative_contrast"
)
audit_schema_ok <- identical(names(denominator_audit), denominator_audit_columns)
audit_rows_ok <- audit_schema_ok &&
  nrow(denominator_audit) == nrow(SCENARIOS) &&
  identical(denominator_audit$scenario_id, SCENARIOS$scenario_id) &&
  identical(denominator_audit$family, SCENARIOS$family) &&
  identical(as.integer(denominator_audit$n), as.integer(SCENARIOS$n)) &&
  close_enough(denominator_audit$rho, SCENARIOS$rho) &&
  all(denominator_audit$paired_replications == CONFIG$n_rep)
audit_numeric_columns <- setdiff(
  denominator_audit_columns, c("scenario_id", "family")
)
audit_finite_ok <- audit_schema_ok && all(vapply(
  denominator_audit[audit_numeric_columns],
  function(column) all(is.finite(column)),
  logical(1)
))
audit_denominators_ok <- audit_finite_ok &&
  all(denominator_audit$nonfinite_cv_test_deviances == 0L) &&
  all(denominator_audit$nonpositive_cv_test_deviances == 0L) &&
  all(denominator_audit$minimum_cv_test_deviance > 0) &&
  all(
    denominator_audit$maximum_cv_test_deviance >=
      denominator_audit$minimum_cv_test_deviance
  ) &&
  all(denominator_audit$nonfinite_relative_contrasts == 0L)
audit_moments_ok <- audit_finite_ok &&
  all(denominator_audit$relative_contrast_variance >= 0) &&
  close_enough(
    denominator_audit$relative_contrast_second_moment,
    denominator_audit$relative_contrast_mean^2 +
      denominator_audit$relative_contrast_variance *
      (denominator_audit$paired_replications - 1) /
      denominator_audit$paired_replications
  ) &&
  all(
    denominator_audit$maximum_absolute_relative_contrast + 5e-13 >=
      abs(denominator_audit$relative_contrast_mean)
  ) &&
  close_enough(
    denominator_audit$relative_contrast_mean,
    scenario$relative_test_deviance_difference
  ) &&
  close_enough(
    sqrt(
      denominator_audit$relative_contrast_variance /
        denominator_audit$paired_replications
    ),
    scenario$relative_test_deviance_difference_mcse
  )
record_check(
  "relative-deviance denominator and moment audit",
  audit_rows_ok && audit_denominators_ok && audit_moments_ok,
  if (audit_finite_ok && nrow(denominator_audit)) {
    paste0(
      "minimum_CV_deviance=",
      format(
        min(denominator_audit$minimum_cv_test_deviance),
        digits = 10L
      ),
      ";maximum_second_moment=",
      format(
        max(denominator_audit$relative_contrast_second_moment),
        scientific = TRUE
      )
    )
  } else {
    "missing, malformed, or nonfinite audit values"
  }
)

margins <- read_result("prediction_margin_sensitivity.csv")
expected_margin_flags <- as.integer(
  prediction_upper < as.numeric(margins$margin)
)
record_check(
  "prediction margin roles and diagnostic flags",
  identical(as.numeric(margins$margin), c(0.001, 0.0005, 0.0025, 0.005)) &&
    identical(margins$role, c("primary", rep("sensitivity", 3L))) &&
    identical(as.integer(margins$noninferiority_pass), expected_margin_flags) &&
    all(abs(margins$upper_one_sided_95 - prediction_upper) < 5e-13),
  paste0("flags=", paste(expected_margin_flags, collapse = ","))
)

diagnostics <- read_result("diagnostic_summary.csv")
zero_diagnostic_columns <- c(
  "full_path_failures", "calibration_failures",
  "calibration_default_uses", "cv_failures",
  "kappa_lower_endpoint_selections", "kappa_upper_endpoint_selections",
  "solver_warnings"
)
diagnostic_required_columns <- c(
  "scenario_id", "attempted_replications", zero_diagnostic_columns,
  "projections", "maximum_radius_interpolation_error"
)
diagnostic_schema_ok <- all(diagnostic_required_columns %in% names(diagnostics))
diagnostic_zero <- diagnostic_schema_ok && all(vapply(
  zero_diagnostic_columns, function(column) {
    all(diagnostics[[column]] == 0L)
  }, logical(1)))
record_check(
  "production numerical diagnostics",
  diagnostic_schema_ok && nrow(diagnostics) == nrow(SCENARIOS) &&
    identical(diagnostics$scenario_id, SCENARIOS$scenario_id) &&
    all(diagnostics$attempted_replications == CONFIG$n_rep) &&
    diagnostic_zero &&
    sum(diagnostics$projections) == 0L &&
    max(diagnostics$maximum_radius_interpolation_error) <
      CONFIG$interpolation_radius_tolerance,
  paste0("projections=", sum(diagnostics$projections),
         ";max_interpolation_error=",
         format(max(diagnostics$maximum_radius_interpolation_error),
                scientific = TRUE))
)

endpoint_methods <- summary$method %in%
  c("PanIC-CF", "CV")
record_check(
  "primary and CV selections avoid radius-grid endpoints",
  all(summary$selected_lower_endpoint_mean[endpoint_methods] == 0) &&
    all(summary$selected_upper_endpoint_mean[endpoint_methods] == 0),
  "PanIC-CF and CV"
)

calibration <- read_result("calibration_target_summary.csv")
record_check(
  "calibration summary is complete",
  nrow(calibration) == nrow(SCENARIOS) &&
    identical(calibration$scenario_id, SCENARIOS$scenario_id) &&
    all(calibration$complete_replications == CONFIG$n_rep),
  paste0("rows=", nrow(calibration))
)

grid <- read_result("grid_sensitivity_summary.csv")
grid_configuration <- read_result("grid_sensitivity_configuration.csv")
grid_config_methods <- grid_configuration$value[
  grid_configuration$item == "methods"
]
grid_required_columns <- c(
  "radius_points", "method", "attempted_replications", "method_failures",
  "full_path_failures", "calibration_failures",
  "calibration_default_uses", "cv_failures",
  "maximum_radius_interpolation_error"
)
grid_schema_ok <- all(grid_required_columns %in% names(grid))
grid_sets_ok <- all(vapply(split(grid$method, grid$radius_points), function(x) {
  setequal(x, assessed_methods)
}, logical(1)))
record_check(
  "grid sensitivity summary and diagnostics",
  grid_schema_ok && nrow(grid) == 9L &&
    identical(sort(unique(grid$radius_points)), c(61L, 121L, 241L)) &&
    grid_sets_ok && all(grid$attempted_replications == 500L) &&
    all(grid$method_failures == 0L) &&
    all(grid$full_path_failures == 0L) &&
    all(grid$calibration_failures == 0L) &&
    all(grid$calibration_default_uses == 0L) &&
    all(grid$cv_failures == 0L) &&
    length(grid_config_methods) == 1L &&
    identical(strsplit(grid_config_methods, ";", fixed = TRUE)[[1L]],
              assessed_methods) &&
    !contains_obsolete_method_label(grid_configuration) &&
    !contains_obsolete_method_label(grid) &&
    max(grid$maximum_radius_interpolation_error) <
      CONFIG$interpolation_radius_tolerance,
  paste0("rows=", nrow(grid), ";max_interpolation_error=",
         format(max(grid$maximum_radius_interpolation_error), scientific = TRUE))
)

grid_pairs <- read_result("grid_sensitivity_paired_contrasts.csv")
grid_contrasts <- read_result("grid_sensitivity_method_contrasts.csv")
panic_grid_contrasts <- grid_contrasts[
  grid_contrasts$lhs_method == "PanIC-CF" &
    grid_contrasts$rhs_method == "CV",
]
grid_pair_keys <- paste(
  grid_pairs$method, grid_pairs$lower_radius_points,
  grid_pairs$upper_radius_points, sep = "\r"
)
grid_point_pairs <- combn(
  sort(CONFIG$grid_sensitivity_points), 2L, simplify = FALSE
)
expected_grid_pair_keys <- unlist(lapply(assessed_methods, function(method) {
  vapply(grid_point_pairs, function(points) {
    paste(method, points[[1L]], points[[2L]], sep = "\r")
  }, character(1))
}), use.names = FALSE)
grid_method_keys <- paste(
  grid_contrasts$radius_points, grid_contrasts$lhs_method,
  grid_contrasts$rhs_method, sep = "\r"
)
expected_grid_method_keys <- unlist(lapply(
  sort(CONFIG$grid_sensitivity_points), function(radius_points) {
    vapply(strsplit(expected_pair_keys, "\r", fixed = TRUE), function(pair) {
      paste(radius_points, pair[[1L]], pair[[2L]], sep = "\r")
    }, character(1))
  }
), use.names = FALSE)
grid_pairs_finite <- all(vapply(
  grid_pairs[vapply(grid_pairs, is.numeric, logical(1))],
  function(column) all(is.finite(column)), logical(1)
))
grid_methods_finite <- all(vapply(
  grid_contrasts[vapply(grid_contrasts, is.numeric, logical(1))],
  function(column) all(is.finite(column)), logical(1)
))
record_check(
  "compact grid contrasts are complete and finite",
  nrow(grid_pairs) == length(expected_grid_pair_keys) &&
    !anyDuplicated(grid_pair_keys) &&
    setequal(grid_pair_keys, expected_grid_pair_keys) &&
    all(
      grid_pairs$paired_replications ==
        CONFIG$grid_sensitivity_replications
    ) &&
    nrow(grid_contrasts) == length(expected_grid_method_keys) &&
    !anyDuplicated(grid_method_keys) &&
    setequal(grid_method_keys, expected_grid_method_keys) &&
    nrow(panic_grid_contrasts) == length(CONFIG$grid_sensitivity_points) &&
    all(
      panic_grid_contrasts$paired_replications ==
        CONFIG$grid_sensitivity_replications
    ) &&
    grid_pairs_finite && grid_methods_finite &&
    !contains_obsolete_method_label(grid_pairs) &&
    !contains_obsolete_method_label(grid_contrasts),
  paste0(
    "grid_pairs=", length(unique(grid_pair_keys)), "/",
    length(expected_grid_pair_keys), ";method_contrasts=",
    length(unique(grid_method_keys)), "/",
    length(expected_grid_method_keys), ";finite=",
    grid_pairs_finite && grid_methods_finite
  )
)

boundary <- read_result("boundary_summary.csv")
subsampling <- read_result("subsampling_summary.csv")
record_check(
  "boundary and subsampling sample sizes match the revised design",
  identical(sort(as.integer(boundary$n)), c(500L, 1000L)) &&
    all(as.integer(subsampling$n) == 1000L) &&
    identical(sort(as.integer(subsampling$b)), c(50L, 100L, 200L)) &&
    close_enough(
      subsampling$b_over_n[order(subsampling$b)], c(0.05, 0.10, 0.20)
    ),
  paste0(
    "boundary_n=", paste(sort(boundary$n), collapse = ","),
    ";subsampling_n=", paste(unique(subsampling$n), collapse = ","),
    ";b=", paste(sort(subsampling$b), collapse = ",")
  )
)

manifest <- read_result("production_file_manifest.csv")
required_manifest_columns <- c("relative_path", "bytes", "sha256")
manifest_paths <- as.character(manifest$relative_path)
manifest_basenames <- basename(manifest_paths)
scenario_artifact_suffixes <- c(
  "_primary.csv", "_diagnostics.csv", "_calibration_rows.csv",
  "_raw.rds", "_runtime.txt"
)
expected_scenario_artifacts <- as.vector(outer(
  SCENARIOS$scenario_id, scenario_artifact_suffixes, paste0
))
scenario_artifact_pattern <- paste0(
  "(", paste(
    gsub("[.]", "[.]", scenario_artifact_suffixes), collapse = "|"
  ), ")$"
)
observed_scenario_artifacts <- manifest_basenames[
  grepl(scenario_artifact_pattern, manifest_basenames) &
    !grepl("^grid_sensitivity_", manifest_basenames)
]
scenario_manifest_ok <-
  !anyDuplicated(observed_scenario_artifacts) &&
  identical(
    sort(observed_scenario_artifacts), sort(expected_scenario_artifacts)
  )
forbidden_manifest_entries <- manifest_basenames[
  grepl("n2000", manifest_basenames, fixed = TRUE) |
    manifest_basenames == "table_confirmatory_decision.tex"
]
record_check(
  "production-run manifest has exactly the ten scenario artifact families",
  identical(names(manifest), required_manifest_columns) &&
    !anyDuplicated(manifest_paths) &&
    all(is.finite(manifest$bytes) & manifest$bytes > 0) &&
    all(grepl("^[0-9a-f]{64}$", manifest$sha256)) &&
    scenario_manifest_ok && !length(forbidden_manifest_entries),
  paste0(
    "entries=", nrow(manifest), ";scenario_artifacts=",
    length(observed_scenario_artifacts), "/",
    length(expected_scenario_artifacts), ";missing=",
    paste(
      setdiff(expected_scenario_artifacts, observed_scenario_artifacts),
      collapse = ","
    ), ";extra=", paste(
      setdiff(observed_scenario_artifacts, expected_scenario_artifacts),
      collapse = ","
    ), ";forbidden=", paste(forbidden_manifest_entries, collapse = ",")
  )
)
active_compact_outputs <- list(
  implementation = implementation,
  configuration = configuration,
  registry = registry,
  seed_ledger = seed_ledger,
  simulation_summary = summary,
  paired_method_contrasts = paired,
  paired_panic_cf_vs_cv = primary_pair,
  scenario_primary_estimands = scenario,
  confirmatory_decision = decision,
  denominator_audit = denominator_audit,
  margin_sensitivity = margins,
  diagnostic_summary = diagnostics,
  calibration_summary = calibration,
  grid_configuration = grid_configuration,
  grid_summary = grid,
  grid_pairs = grid_pairs,
  grid_method_contrasts = grid_contrasts,
  boundary_summary = boundary,
  subsampling_summary = subsampling,
  production_manifest = manifest
)
obsolete_output_frames <- names(active_compact_outputs)[vapply(
  active_compact_outputs, contains_obsolete_method_label, logical(1)
)]
record_check(
  "active compact outputs contain no obsolete method labels",
  !length(obsolete_output_frames),
  if (length(obsolete_output_frames)) {
    paste(obsolete_output_frames, collapse = ",")
  } else {
    "none"
  }
)

table_names <- PANIC_RENDERED_TABLE_FILES
verify_table_mirrors <- function() {
  temporary_dir <- tempfile(pattern = "panic-rendered-tables-")
  dir.create(temporary_dir)
  on.exit(unlink(temporary_dir, recursive = TRUE, force = TRUE), add = TRUE)
  rendered <- render_all_manuscript_tables(results_dir)
  write_all_manuscript_tables(rendered, temporary_dir)
  matches <- vapply(table_names, function(name) {
    expected_path <- file.path(results_dir, name)
    rendered_path <- file.path(temporary_dir, name)
    if (!file.exists(expected_path) || !file.exists(rendered_path)) {
      return(FALSE)
    }
    expected_size <- file.info(expected_path)$size
    rendered_size <- file.info(rendered_path)$size
    identical(expected_size, rendered_size) &&
      identical(
        readBin(expected_path, what = "raw", n = expected_size),
        readBin(rendered_path, what = "raw", n = rendered_size)
      )
  }, logical(1))
  list(matches = matches, error = NULL)
}
table_mirror_result <- tryCatch(
  verify_table_mirrors(),
  error = function(error) {
    list(
      matches = setNames(rep(FALSE, length(table_names)), table_names),
      error = conditionMessage(error)
    )
  }
)
record_check(
  "all six canonical table mirrors reproduce byte-for-byte",
  all(table_mirror_result$matches),
  if (is.null(table_mirror_result$error)) {
    paste0(
      "matched=", sum(table_mirror_result$matches), "/", length(table_names)
    )
  } else {
    table_mirror_result$error
  }
)
manuscript_path <- file.path(repository_dir, "manuscript", "main.tex")
inline_table_result <- tryCatch(
  verify_inline_tables(
    manuscript_path, results_dir, stop_on_failure = FALSE
  ),
  error = function(error) error
)
inline_tables_ok <- !inherits(inline_table_result, "error") &&
  all(inline_table_result$passed == 1L) &&
  isTRUE(attr(inline_table_result, "no_external_table_references")) &&
  isTRUE(attr(inline_table_result, "exact_marker_inventory"))
record_check(
  "manuscript embeds all six mirrors with no external table dependency",
  inline_tables_ok,
  if (inherits(inline_table_result, "error")) {
    conditionMessage(inline_table_result)
  } else {
    paste0(
      "matched=", sum(inline_table_result$passed), "/",
      nrow(inline_table_result), ";external_references=0;marker_inventory=",
      if (isTRUE(attr(inline_table_result, "exact_marker_inventory"))) {
        "exact"
      } else {
        attr(inline_table_result, "marker_inventory_detail")
      }
    )
  }
)
table_text <- paste(vapply(table_names, function(name) {
  paste(readLines(assert_file(name), warn = FALSE), collapse = "\n")
}, character(1)), collapse = "\n")
record_check(
  "manuscript-facing tables use the three assessed labels",
  !grepl(
    paste0(
      "TODO|", obsolete_method_pattern,
      "|second-confirmation|revised PanIC-CF"
    ),
    table_text,
    ignore.case = TRUE
  ) &&
    grepl(" & PanIC-CF & ", table_text, fixed = TRUE) &&
    grepl(" & BIC-like & ", table_text, fixed = TRUE) &&
    grepl(" & CV & ", table_text, fixed = TRUE),
  paste(table_names, collapse = ",")
)
manuscript_text <- paste(readLines(manuscript_path, warn = FALSE),
                         collapse = "\n")
record_check(
  "manuscript contains no TODO or obsolete assessed-method labels",
  !grepl("TODO", manuscript_text, ignore.case = TRUE) &&
    !grepl(obsolete_method_pattern, manuscript_text),
  "manuscript/main.tex"
)

artifact_names <- c("figure_runtime.pdf", "figure_boundary_subsampling.pdf")
artifact_paths <- file.path(results_dir, artifact_names)
record_check(
  "manuscript figures are present and nonempty",
  all(file.exists(artifact_paths)) && all(file.info(artifact_paths)$size > 1000),
  paste0(artifact_names, "=", file.info(artifact_paths)$size, collapse = ";")
)

verification <- do.call(rbind, checks)
print(verification, row.names = FALSE)
if (any(verification$passed != 1L)) {
  failed <- verification$check[verification$passed != 1L]
  stop("Repository verification failed: ", paste(failed, collapse = "; "),
       call. = FALSE)
}
write.csv(
  verification,
  file.path(results_dir, "repository_verification.csv"),
  row.names = FALSE
)
cat("PASS: compact PanIC-CF repository release is verified.\n")
