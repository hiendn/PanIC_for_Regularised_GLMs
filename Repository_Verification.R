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
    output <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  } else if (nzchar(Sys.which("sha256sum"))) {
    output <- system2("sha256sum", path, stdout = TRUE)
  } else {
    stop("Neither shasum nor sha256sum is available", call. = FALSE)
  }
  substr(output[[1L]], 1L, 64L)
}
close_enough <- function(x, y, tolerance = 5e-13) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance,
                   check.attributes = FALSE))
}

## Check the compact release checksum manifest first.
checksum_path <- assert_file("RESULTS_SHA256.txt")
checksum_lines <- readLines(checksum_path, warn = FALSE)
checksum_entries <- grep("^[0-9a-f]{64}  ", checksum_lines, value = TRUE)
checksum_expected <- substr(checksum_entries, 1L, 64L)
checksum_relative <- substring(checksum_entries, 67L)
checksum_files <- file.path(repository_dir, checksum_relative)
checksum_present <- length(checksum_entries) > 0L && all(file.exists(checksum_files))
checksum_actual <- if (checksum_present) {
  vapply(checksum_files, sha256_file, character(1))
} else {
  rep(NA_character_, length(checksum_expected))
}
record_check(
  "compact result checksum manifest",
  checksum_present && identical(unname(checksum_actual), checksum_expected),
  paste0("entries=", length(checksum_entries))
)

implementation <- read_result("implementation_validation.csv")
record_check(
  "all deterministic implementation checks passed",
  nrow(implementation) == 21L && all(implementation$passed == 1L),
  paste0("passed=", sum(implementation$passed), "/", nrow(implementation))
)

configuration <- read_result("configuration.csv")
config_value <- function(item) {
  value <- configuration$value[configuration$item == item]
  if (length(value) != 1L) stop("Configuration item is not unique: ", item)
  value
}
record_check(
  "production configuration and decision constants",
  as.integer(config_value("master_seed")) == 2066091802L &&
    as.integer(config_value("replications_requested")) == 1000L &&
    as.integer(config_value("independent_test_size")) == 2000L &&
    as.numeric(config_value("prediction_noninferiority_margin")) == 0.001 &&
    config_value("cv_primary_rule") ==
      "minimum mean five-fold validation loss on radius grid",
  paste0("seed=", config_value("master_seed"),
         ";replications=", config_value("replications_requested"))
)

registry <- read_result("scenario_registry.csv")
record_check(
  "seven locked scenarios are present",
  nrow(registry) == 7L &&
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
seed_values_ok <- seed_schema_ok && nrow(seed_ledger) == 7000L
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
expected_methods <- c(
  "PanIC-CF", "PanIC-CF-original", "BIC-like",
  "BIC-active (exploratory)", "CV-min", "CV-1SE"
)
method_sets_ok <- all(vapply(split(summary, summary$scenario_id), function(dat) {
  expected <- c(
    "PanIC-CF", "PanIC-CF-original",
    if (dat$family[1L] == "gaussian") "BIC-like" else
      "BIC-active (exploratory)",
    "CV-min", "CV-1SE"
  )
  setequal(dat$method, expected)
}, logical(1)))
record_check(
  "primary summary is complete",
  nrow(summary) == 35L && method_sets_ok &&
    all(summary$attempted_replications == 1000L) &&
    all(summary$failed_replications == 0L) &&
    all(summary$method %in% expected_methods),
  paste0("rows=", nrow(summary), ";failures=",
         sum(summary$failed_replications))
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
decision_values_ok <- nrow(scenario) == 7L && nrow(decision) == 1L &&
  all(scenario$paired_replications == 1000L) &&
  close_enough(decision$support_estimate, support_estimate) &&
  close_enough(decision$support_mcse, support_mcse) &&
  close_enough(decision$support_upper_one_sided_95, support_upper) &&
  close_enough(decision$prediction_estimate, prediction_estimate) &&
  close_enough(decision$prediction_mcse, prediction_mcse) &&
  close_enough(decision$prediction_upper_one_sided_95, prediction_upper)
record_check(
  "confirmatory endpoints reconstruct from scenario rows",
  decision_values_ok,
  sprintf("support=%.9f;relative_deviance=%.9f",
          support_estimate, prediction_estimate)
)
record_check(
  "locked joint superiority/noninferiority decision passes",
  decision_values_ok && support_upper < 0 && prediction_upper < 0.001 &&
    decision$complete_pairing == 1L &&
    decision$support_superiority_pass == 1L &&
    decision$prediction_noninferiority_pass == 1L &&
    decision$joint_claim_pass == 1L,
  sprintf("support_upper=%.9f;prediction_upper=%.9f",
          support_upper, prediction_upper)
)

denominator_audit <- read_result("relative_deviance_denominator_audit.csv")
denominator_audit_columns <- c(
  "scenario_id", "family", "n", "rho", "paired_replications",
  "nonfinite_cv_min_test_deviances",
  "nonpositive_cv_min_test_deviances",
  "minimum_cv_min_test_deviance",
  "maximum_cv_min_test_deviance",
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
  all(denominator_audit$paired_replications == 1000L)
audit_numeric_columns <- setdiff(
  denominator_audit_columns, c("scenario_id", "family")
)
audit_finite_ok <- audit_schema_ok && all(vapply(
  denominator_audit[audit_numeric_columns],
  function(column) all(is.finite(column)),
  logical(1)
))
audit_denominators_ok <- audit_finite_ok &&
  all(denominator_audit$nonfinite_cv_min_test_deviances == 0L) &&
  all(denominator_audit$nonpositive_cv_min_test_deviances == 0L) &&
  all(denominator_audit$minimum_cv_min_test_deviance > 0) &&
  all(
    denominator_audit$maximum_cv_min_test_deviance >=
      denominator_audit$minimum_cv_min_test_deviance
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
      "minimum_CV_min_deviance=",
      format(
        min(denominator_audit$minimum_cv_min_test_deviance),
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
record_check(
  "prediction margin roles and decisions",
  identical(as.numeric(margins$margin), c(0.001, 0.0005, 0.0025, 0.005)) &&
    identical(margins$role, c("primary", rep("sensitivity", 3L))) &&
    all(margins$noninferiority_pass == 1L) &&
    all(abs(margins$upper_one_sided_95 - prediction_upper) < 5e-13),
  paste0("smallest_passing_margin=", min(
    margins$margin[margins$noninferiority_pass == 1L]
  ))
)

diagnostics <- read_result("diagnostic_summary.csv")
zero_diagnostic_columns <- c(
  "full_path_failures", "calibration_failures",
  "calibration_default_uses", "cv_failures",
  "kappa_lower_endpoint_selections", "kappa_upper_endpoint_selections",
  "original_kappa_lower_endpoint_selections",
  "original_kappa_upper_endpoint_selections", "solver_warnings"
)
diagnostic_zero <- all(vapply(zero_diagnostic_columns, function(column) {
  all(diagnostics[[column]] == 0L)
}, logical(1)))
record_check(
  "production numerical diagnostics",
  nrow(diagnostics) == 7L &&
    all(diagnostics$attempted_replications == 1000L) && diagnostic_zero &&
    sum(diagnostics$projections) == 1L &&
    max(diagnostics$maximum_radius_interpolation_error) <
      CONFIG$interpolation_radius_tolerance,
  paste0("projections=", sum(diagnostics$projections),
         ";max_interpolation_error=",
         format(max(diagnostics$maximum_radius_interpolation_error),
                scientific = TRUE))
)

endpoint_methods <- summary$method %in%
  c("PanIC-CF", "PanIC-CF-original", "CV-min", "CV-1SE")
record_check(
  "primary and CV selections avoid radius-grid endpoints",
  all(summary$selected_lower_endpoint_mean[endpoint_methods] == 0) &&
    all(summary$selected_upper_endpoint_mean[endpoint_methods] == 0),
  "PanIC-CF, original-weight PanIC-CF, CV-min, and CV-1SE"
)

calibration <- read_result("calibration_target_summary.csv")
record_check(
  "calibration summary is complete",
  nrow(calibration) == 7L &&
    identical(calibration$scenario_id, SCENARIOS$scenario_id) &&
    all(calibration$complete_replications == 1000L),
  paste0("rows=", nrow(calibration))
)

grid <- read_result("grid_sensitivity_summary.csv")
grid_expected_methods <- c(
  "PanIC-CF", "PanIC-CF-original", "BIC-like", "CV-min", "CV-1SE"
)
grid_sets_ok <- all(vapply(split(grid$method, grid$radius_points), function(x) {
  setequal(x, grid_expected_methods)
}, logical(1)))
record_check(
  "grid sensitivity summary and diagnostics",
  nrow(grid) == 15L &&
    identical(sort(unique(grid$radius_points)), c(61L, 121L, 241L)) &&
    grid_sets_ok && all(grid$attempted_replications == 500L) &&
    all(grid$method_failures == 0L) &&
    all(grid$full_path_failures == 0L) &&
    all(grid$calibration_failures == 0L) &&
    all(grid$calibration_default_uses == 0L) &&
    all(grid$cv_failures == 0L) &&
    max(grid$maximum_radius_interpolation_error) <
      CONFIG$interpolation_radius_tolerance,
  paste0("rows=", nrow(grid), ";max_interpolation_error=",
         format(max(grid$maximum_radius_interpolation_error), scientific = TRUE))
)

grid_contrasts <- read_result("grid_sensitivity_method_contrasts.csv")
panic_grid_contrasts <- grid_contrasts[
  grid_contrasts$lhs_method == "PanIC-CF" &
    grid_contrasts$rhs_method == "CV-min",
]
record_check(
  "PanIC-CF grid contrasts to CV-min are present and favorable",
  nrow(panic_grid_contrasts) == 3L &&
    all(panic_grid_contrasts$paired_replications == 500L) &&
    all(panic_grid_contrasts$wrong_difference < 0),
  paste(sprintf("m=%d:%.3f", panic_grid_contrasts$radius_points,
                panic_grid_contrasts$wrong_difference), collapse = ";")
)

manifest <- read_result("production_file_manifest.csv")
required_manifest_columns <- c("relative_path", "bytes", "sha256")
record_check(
  "historical production-run file manifest is well formed",
  all(required_manifest_columns %in% names(manifest)) && nrow(manifest) >= 70L &&
    all(nchar(manifest$sha256) == 64L) &&
    any(grepl("_raw[.]rds$", manifest$relative_path)) &&
    any(grepl("_calibration_rows[.]csv$", manifest$relative_path)),
  paste0("entries=", nrow(manifest))
)

table_names <- c(
  "table_primary_support.tex", "table_primary_performance.tex",
  "table_confirmatory_decision.tex", "table_calibration.tex",
  "table_grid_sensitivity.tex", "table_boundary.tex", "table_subsampling.tex"
)
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
  "all seven canonical table mirrors reproduce byte-for-byte",
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
  isTRUE(attr(inline_table_result, "no_external_table_references"))
record_check(
  "manuscript embeds all seven mirrors with no external table dependency",
  inline_tables_ok,
  if (inherits(inline_table_result, "error")) {
    conditionMessage(inline_table_result)
  } else {
    paste0(
      "matched=", sum(inline_table_result$passed), "/",
      nrow(inline_table_result), ";external_references=0"
    )
  }
)
table_text <- paste(vapply(table_names, function(name) {
  paste(readLines(assert_file(name), warn = FALSE), collapse = "\n")
}, character(1)), collapse = "\n")
decision_table <- paste(
  readLines(assert_file("table_confirmatory_decision.tex"), warn = FALSE),
  collapse = "\n"
)
record_check(
  "manuscript-facing tables use final labels and no TODO markers",
  !grepl("TODO|BIC-active|second-confirmation|revised PanIC-CF", table_text,
         ignore.case = TRUE) &&
    grepl("CV-min", table_text, fixed = TRUE) &&
    grepl("CV-1SE", table_text, fixed = TRUE) &&
    grepl("exploratory", table_text, ignore.case = TRUE),
  paste(table_names, collapse = ",")
)
record_check(
  "confirmatory table agrees with locked decision",
  grepl("-0.43814", decision_table, fixed = TRUE) &&
    grepl("-0.40854", decision_table, fixed = TRUE) &&
    grepl("0.000160", decision_table, fixed = TRUE) &&
    grepl("0.000214", decision_table, fixed = TRUE) &&
    grepl("<0.001", gsub("\\\\", "", decision_table), fixed = TRUE),
  "displayed estimates, upper bounds, and decision thresholds"
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
