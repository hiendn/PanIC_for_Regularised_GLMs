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
  normalizePath("Results_Verification.R")
}
candidate_dir <- dirname(script_path)
source(file.path(candidate_dir, "Method_Lock_Verification.R"))
verify_method_lock(candidate_dir)
source(file.path(candidate_dir, "Simulation_Config.R"))
source(file.path(candidate_dir, "Manuscript_Table_Rendering.R"))
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
  assert(file.exists(path), paste("Missing result:", name))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
mcse <- function(x) sd(x) / sqrt(length(x))

configuration <- read_result("configuration.csv")
n_rep <- as.integer(configuration$value[
  configuration$item == "replications_requested"
])
assert(length(n_rep) == 1L && is.finite(n_rep),
       "Cannot recover the requested replication count")
seed_role <- configuration$value[configuration$item == "seed_role"]
assert(length(seed_role) == 1L &&
         seed_role %in% c("production", "dedicated-smoke"),
       "The seed role is missing or invalid")
expected_master_seed <- if (seed_role == "production") {
  CONFIG$master_seed
} else {
  CONFIG$smoke_master_seed
}
assert(as.integer(configuration$value[
  configuration$item == "master_seed"
]) == expected_master_seed, "The result master seed is not locked")
assert(as.integer(configuration$value[
  configuration$item == "production_master_seed"
]) == CONFIG$master_seed, "The recorded production seed has changed")
if (seed_role == "production") {
  assert(n_rep == CONFIG$n_rep,
         "The production archive does not contain the locked replication count")
}

validation <- read.csv(
  file.path(candidate_dir, "results", "implementation_validation.csv"),
                       stringsAsFactors = FALSE)
assert(all(validation$passed == 1L),
       "A deterministic implementation check failed")

expected_methods <- function(family) {
  c("PanIC-CF", "BIC-like", CONFIG$primary_cv_method)
}
all_primary <- list()
all_diagnostics <- list()
all_calibration <- list()
for (scenario_index in seq_len(nrow(SCENARIOS))) {
  id <- SCENARIOS$scenario_id[scenario_index]
  primary <- read_result(paste0(id, "_primary.csv"))
  diagnostic <- read_result(paste0(id, "_diagnostics.csv"))
  calibration <- read_result(paste0(id, "_calibration_rows.csv"))
  family <- SCENARIOS$family[scenario_index]
  all_primary[[id]] <- primary
  all_diagnostics[[id]] <- diagnostic
  all_calibration[[id]] <- calibration

  assert(nrow(primary) == 3L * n_rep,
         paste0(id, ": expected three method rows per replication"))
  assert(nrow(diagnostic) == n_rep,
         paste0(id, ": diagnostic row count mismatch"))
  assert(nrow(calibration) == 10L * n_rep,
         paste0(id, ": calibration row count mismatch"))
  assert(identical(sort(unique(diagnostic$replication)), seq_len(n_rep)),
         paste0(id, ": incomplete replication sequence"))
  assert(all(table(primary$replication) == 3L),
         paste0(id, ": a replication does not contain three methods"))
  assert(all(vapply(
    split(primary$method, primary$replication),
    function(methods) setequal(methods, expected_methods(family)),
    logical(1)
  )), paste0(id, ": method set mismatch"))
  assert(all(table(calibration$replication) == 10L),
         paste0(id, ": calibration rows are incomplete"))

  ok_target <- is.finite(calibration$raw_cross_signed_radius)
  if (any(ok_target)) {
    expected_primary_weight <- sqrt(log(log(
      calibration$n_validation[ok_target] + exp(exp(1))
    )))
    assert(max(abs(
      expected_primary_weight - calibration$weight[ok_target]
    )) < 1e-12, paste0(id, ": calibration weight mismatch"))
  }

  pan <- primary[primary$method == "PanIC-CF", ]
  cv <- primary[primary$method == CONFIG$primary_cv_method, ]
  assert(all(pan$calibration_failed == diagnostic$calibration_failed),
         paste0(id, ": calibration failure flag mismatch"))
  assert(all(pan$default_used == diagnostic$default_used),
         paste0(id, ": calibration default flag mismatch"))
  assert(all(cv$method_failed == diagnostic$cv_failed),
         paste0(id, ": CV failure flag mismatch"))

  good_path <- diagnostic$full_path_failed == 0L
  for (column in "kappa_hat") {
    on_grid <- vapply(diagnostic[[column]][good_path], function(value) {
      min(abs(value - CONFIG$kappa_grid)) < 1e-12
    }, logical(1))
    assert(all(on_grid), paste0(id, ": ", column, " is off-grid"))
  }
  finite_error <- diagnostic$maximum_active_radius_error[
    is.finite(diagnostic$maximum_active_radius_error)
  ]
  assert(!length(finite_error) ||
           max(finite_error) <= CONFIG$interpolation_radius_tolerance,
         paste0(id, ": interpolation tolerance exceeded"))
}

primary <- do.call(rbind, all_primary)
diagnostics <- do.call(rbind, all_diagnostics)
seed_ledger <- read_result("seed_ledger.csv")
assert(nrow(seed_ledger) == nrow(SCENARIOS) * n_rep,
       "Seed ledger does not cover all attempted replications")
assert(!anyDuplicated(seed_ledger[c("scenario_id", "replication")]),
       "Seed ledger duplicates a scenario/replication pair")
seed_columns <- setdiff(names(seed_ledger), c("scenario_id", "replication"))
assert(!anyDuplicated(unlist(seed_ledger[seed_columns], use.names = FALSE)),
       "A seed is reused anywhere in the main production ledger")
assert(all(diagnostics$calibration_rows_expected == 10L),
       "A diagnostic row reports a non-locked calibration-row count")

simulation_summary <- read_result("simulation_summary.csv")
paired <- read_result("paired_method_contrasts.csv")
scenario_estimands <- read_result("scenario_primary_estimands.csv")
decision <- read_result("confirmatory_decision.csv")
sensitivity <- read_result("prediction_margin_sensitivity.csv")
assert(nrow(simulation_summary) == 3L * nrow(SCENARIOS),
       "Simulation summary must contain three methods in every scenario")
assert(nrow(paired) == 3L * nrow(SCENARIOS),
       "Paired summary must contain three comparisons in every scenario")
assert(nrow(scenario_estimands) == nrow(SCENARIOS) && nrow(decision) == 1L,
       "Scenario estimand outputs have the wrong dimensions")
assert(nrow(sensitivity) == 4L && sum(sensitivity$role == "primary") == 1L,
       "Prediction-margin sensitivity output is incomplete")

## Reconstruct every field used by the primary support and performance tables
## from the replication-level rows.  The later byte check then verifies
## presentation fidelity to a summary that has itself been independently
## reproduced from raw output.
summary_metrics <- c(
  "fp", "fn", "fpr", "fnr", "exact", "wrong",
  "signed_attained_radius_error", "test_deviance",
  "selected_grid_radius", "attained_radius",
  "selected_lower_endpoint", "selected_upper_endpoint"
)
summary_mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x) / sqrt(length(x))
}
expected_summary_rows <- lapply(
  split(primary, interaction(primary$scenario_id, primary$method, drop = TRUE)),
  function(dat) {
    out <- dat[1L, c("scenario_id", "family", "n", "rho", "method")]
    out$attempted_replications <- nrow(dat)
    out$failed_replications <- sum(dat$method_failed)
    out$successful_replications <- sum(dat$method_failed == 0L)
    for (metric in summary_metrics) {
      out[[paste0(metric, "_mean")]] <- mean(dat[[metric]], na.rm = TRUE)
      out[[paste0(metric, "_mcse")]] <- summary_mcse(dat[[metric]])
    }
    out
  }
)
expected_simulation_summary <- do.call(rbind, expected_summary_rows)
reported_key <- paste(
  simulation_summary$scenario_id, simulation_summary$method, sep = "\r"
)
expected_key <- paste(
  expected_simulation_summary$scenario_id,
  expected_simulation_summary$method,
  sep = "\r"
)
summary_match <- match(reported_key, expected_key)
assert(
  !anyNA(summary_match) && !anyDuplicated(reported_key) &&
    !anyDuplicated(expected_key),
  "Simulation-summary scenario/method keys are incomplete or duplicated"
)
expected_simulation_summary <-
  expected_simulation_summary[summary_match, , drop = FALSE]
rownames(expected_simulation_summary) <- NULL
rownames(simulation_summary) <- NULL
summary_character_columns <- c("scenario_id", "family", "method")
summary_numeric_columns <- setdiff(
  names(simulation_summary), summary_character_columns
)
assert(
  identical(
    simulation_summary[summary_character_columns],
    expected_simulation_summary[summary_character_columns]
  ) &&
    isTRUE(all.equal(
      simulation_summary[summary_numeric_columns],
      expected_simulation_summary[summary_numeric_columns],
      tolerance = 1e-12, check.attributes = FALSE
    )),
  "Primary support/performance summary was not reproduced from raw rows"
)

## Independently recompute the locked primary estimands and bounds from raw
## replication-level rows, rather than trusting the analysis output.
scenario_support <- scenario_support_se <- numeric(nrow(SCENARIOS))
scenario_prediction <- scenario_prediction_se <- numeric(nrow(SCENARIOS))
paired_counts <- integer(nrow(SCENARIOS))
denominator_audit_rows <- vector("list", nrow(SCENARIOS))
for (i in seq_len(nrow(SCENARIOS))) {
  dat <- all_primary[[SCENARIOS$scenario_id[i]]]
  lhs <- dat[dat$method == "PanIC-CF" & dat$method_failed == 0L, ]
  rhs <- dat[
    dat$method == CONFIG$primary_cv_method & dat$method_failed == 0L,
  ]
  pair <- merge(lhs, rhs, by = "replication", suffixes = c("_lhs", "_rhs"))
  paired_counts[i] <- nrow(pair)
  support_delta <- pair$wrong_lhs - pair$wrong_rhs
  denominator <- pair$test_deviance_rhs
  assert(
    all(is.finite(denominator) & denominator > 0),
    paste0(
      SCENARIOS$scenario_id[i],
      ": relative test-deviance denominators are not all finite and positive"
    )
  )
  numerator <- pair$test_deviance_lhs - denominator
  assert(
    all(is.finite(numerator)),
    paste0(
      SCENARIOS$scenario_id[i],
      ": relative test-deviance numerators are not all finite"
    )
  )
  prediction_delta <-
    numerator / denominator
  assert(
    all(is.finite(prediction_delta)),
    paste0(
      SCENARIOS$scenario_id[i],
      ": relative test-deviance contrasts are not all finite"
    )
  )
  denominator_audit_rows[[i]] <- data.frame(
    scenario_id = SCENARIOS$scenario_id[i],
    family = SCENARIOS$family[i],
    n = SCENARIOS$n[i],
    rho = SCENARIOS$rho[i],
    paired_replications = nrow(pair),
    nonfinite_cv_test_deviances = sum(!is.finite(denominator)),
    nonpositive_cv_test_deviances = sum(
      is.finite(denominator) & denominator <= 0
    ),
    minimum_cv_test_deviance = min(denominator),
    maximum_cv_test_deviance = max(denominator),
    nonfinite_relative_contrasts = sum(!is.finite(prediction_delta)),
    relative_contrast_mean = mean(prediction_delta),
    relative_contrast_variance = var(prediction_delta),
    relative_contrast_second_moment = mean(prediction_delta^2),
    maximum_absolute_relative_contrast = max(abs(prediction_delta)),
    stringsAsFactors = FALSE
  )
  scenario_support[i] <- mean(support_delta)
  scenario_support_se[i] <- mcse(support_delta)
  scenario_prediction[i] <- mean(prediction_delta)
  scenario_prediction_se[i] <- mcse(prediction_delta)
}
expected_denominator_audit <- do.call(rbind, denominator_audit_rows)
reported_denominator_audit <- read_result(
  "relative_deviance_denominator_audit.csv"
)
assert(
  identical(names(reported_denominator_audit), names(expected_denominator_audit)) &&
    identical(
      reported_denominator_audit$scenario_id,
      expected_denominator_audit$scenario_id
    ) &&
    isTRUE(all.equal(
      reported_denominator_audit[setdiff(
        names(reported_denominator_audit), c("scenario_id", "family")
      )],
      expected_denominator_audit[setdiff(
        names(expected_denominator_audit), c("scenario_id", "family")
      )],
      tolerance = 1e-12, check.attributes = FALSE
    )) &&
    identical(
      reported_denominator_audit$family,
      expected_denominator_audit$family
    ),
  "The relative-deviance denominator audit was not reproduced from raw rows"
)
support_estimate <- mean(scenario_support)
support_mcse <- sqrt(sum(scenario_support_se^2)) / nrow(SCENARIOS)
support_upper <- support_estimate + qnorm(0.95) * support_mcse
prediction_estimate <- mean(scenario_prediction)
prediction_mcse <- sqrt(sum(scenario_prediction_se^2)) / nrow(SCENARIOS)
prediction_upper <- prediction_estimate + qnorm(0.95) * prediction_mcse
close <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-12))
assert(close(decision$support_estimate, support_estimate) &&
         close(decision$support_mcse, support_mcse) &&
         close(decision$support_upper_one_sided_95, support_upper),
       "Support pooled summary was not reproduced from raw rows")
assert(close(decision$prediction_estimate, prediction_estimate) &&
         close(decision$prediction_mcse, prediction_mcse) &&
         close(decision$prediction_upper_one_sided_95, prediction_upper),
       "Prediction pooled summary was not reproduced from raw rows")
expected_complete <- all(paired_counts == n_rep)
expected_support_pass <- is.finite(support_upper) && support_upper < 0
expected_prediction_pass <- is.finite(prediction_upper) &&
  prediction_upper < CONFIG$prediction_noninferiority_margin
assert(decision$complete_pairing == as.integer(expected_complete) &&
         decision$support_superiority_pass ==
           as.integer(expected_support_pass) &&
         decision$prediction_noninferiority_pass ==
           as.integer(expected_prediction_pass) &&
         decision$joint_claim_pass == as.integer(
           expected_complete && expected_support_pass &&
             expected_prediction_pass
         ), "The reported endpoint flags were computed incorrectly")

required_artifacts <- c(
  "diagnostic_summary.csv", "calibration_target_summary.csv",
  "relative_deviance_denominator_audit.csv",
  "table_primary_support.tex", "table_primary_performance.tex",
  "table_calibration.tex",
  "sessionInfo.txt", "environment.txt", "locked_configuration.rds"
)
assert(all(file.exists(file.path(results_dir, required_artifacts))),
       "One or more required simulation artifacts are missing")

reported_diagnostic_summary <- read_result("diagnostic_summary.csv")
reported_calibration_summary <- read_result("calibration_target_summary.csv")
expected_tables <- setNames(
  list(
    render_primary_support_table(simulation_summary, configuration),
    render_primary_performance_table(simulation_summary, configuration),
    render_calibration_table(
      reported_calibration_summary, reported_diagnostic_summary
    )
  ),
  c(
    "table_primary_support.tex", "table_primary_performance.tex",
    "table_calibration.tex"
  )
)
for (table_name in names(expected_tables)) {
  path <- file.path(results_dir, table_name)
  expected_bytes <- panic_table_bytes(expected_tables[[table_name]])
  observed_size <- file.info(path)$size
  observed_bytes <- readBin(path, what = "raw", n = observed_size)
  assert(
    identical(observed_size, as.numeric(length(expected_bytes))) &&
      identical(observed_bytes, expected_bytes),
    paste0(
      "Canonical table mirror is stale or non-deterministic: ", table_name
    )
  )
}
support_table_text <- paste(readLines(
  file.path(results_dir, "table_primary_support.tex"),
  warn = FALSE
), collapse = "\n")
assert(!grepl("BIC-active \\(exploratory\\)", support_table_text) &&
         grepl("BIC-like", support_table_text),
       "Reader-facing BIC-like labels were not normalized")
assert(!grepl(
  "second-confirmation|revised PanIC-CF", support_table_text,
  ignore.case = TRUE
), "A manuscript table exposes internal process wording")

cat("PASS: simulation results and reported endpoint flags are verified.\n")
cat("Pooled diagnostic flag:", decision$joint_claim_pass, "\n")
cat("Support upper bound:", format(support_upper, digits = 8), "\n")
cat("Prediction upper bound:", format(prediction_upper, digits = 8), "\n")
