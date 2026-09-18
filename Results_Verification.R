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
  c(
    "PanIC-CF", CONFIG$original_sensitivity_method,
    if (family == "gaussian") "BIC-like" else "BIC-active (exploratory)",
    CONFIG$primary_cv_method, CONFIG$secondary_cv_method
  )
}
all_primary <- list()
all_diagnostics <- list()
for (scenario_index in seq_len(nrow(SCENARIOS))) {
  id <- SCENARIOS$scenario_id[scenario_index]
  primary <- read_result(paste0(id, "_primary.csv"))
  diagnostic <- read_result(paste0(id, "_diagnostics.csv"))
  calibration <- read_result(paste0(id, "_calibration_rows.csv"))
  family <- SCENARIOS$family[scenario_index]
  all_primary[[id]] <- primary
  all_diagnostics[[id]] <- diagnostic

  assert(nrow(primary) == 5L * n_rep,
         paste0(id, ": expected five method rows per replication"))
  assert(nrow(diagnostic) == n_rep,
         paste0(id, ": diagnostic row count mismatch"))
  assert(nrow(calibration) == 10L * n_rep,
         paste0(id, ": calibration row count mismatch"))
  assert(identical(sort(unique(diagnostic$replication)), seq_len(n_rep)),
         paste0(id, ": incomplete replication sequence"))
  assert(all(table(primary$replication) == 5L),
         paste0(id, ": a replication does not contain five methods"))
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
    expected_original_weight <- log(log(
      calibration$n_validation[ok_target] + exp(exp(1))
    ))
    assert(max(abs(
      expected_primary_weight - calibration$weight[ok_target]
    )) < 1e-12, paste0(id, ": revised weight mismatch"))
    assert(max(abs(
      expected_original_weight - calibration$original_weight[ok_target]
    )) < 1e-12, paste0(id, ": original sensitivity weight mismatch"))
  }

  pan <- primary[primary$method == "PanIC-CF", ]
  original <- primary[
    primary$method == CONFIG$original_sensitivity_method,
  ]
  cv_min <- primary[primary$method == CONFIG$primary_cv_method, ]
  cv_one_se <- primary[primary$method == CONFIG$secondary_cv_method, ]
  assert(all(pan$calibration_failed == diagnostic$calibration_failed) &&
           all(original$calibration_failed == diagnostic$calibration_failed),
         paste0(id, ": calibration failure flag mismatch"))
  assert(all(pan$default_used == diagnostic$default_used) &&
           all(original$default_used == diagnostic$default_used),
         paste0(id, ": calibration default flag mismatch"))
  assert(all(cv_min$method_failed == diagnostic$cv_failed) &&
           all(cv_one_se$method_failed == diagnostic$cv_failed),
         paste0(id, ": CV failure flag mismatch"))
  cv_ok <- cv_min$method_failed == 0L & cv_one_se$method_failed == 0L
  assert(all(
    cv_one_se$selected_grid_radius[cv_ok] <=
      cv_min$selected_grid_radius[cv_ok] + 1e-14
  ), paste0(id, ": CV-1SE selected a radius above CV-min"))
  assert(all(
    diagnostic$cv_one_se_index[diagnostic$cv_failed == 0L] <=
      diagnostic$cv_min_index[diagnostic$cv_failed == 0L]
  ), paste0(id, ": CV index ordering mismatch"))

  good_path <- diagnostic$full_path_failed == 0L
  for (column in c("kappa_hat", "original_kappa_hat")) {
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
       "A seed is reused anywhere in the main confirmatory ledger")
assert(all(diagnostics$calibration_rows_expected == 10L),
       "A diagnostic row reports a non-locked calibration-row count")

simulation_summary <- read_result("simulation_summary.csv")
paired <- read_result("paired_method_contrasts.csv")
scenario_estimands <- read_result("scenario_primary_estimands.csv")
decision <- read_result("confirmatory_decision.csv")
sensitivity <- read_result("prediction_margin_sensitivity.csv")
assert(nrow(simulation_summary) == 35L,
       "Simulation summary must contain five methods in seven scenarios")
assert(nrow(paired) == 42L,
       "Paired summary must contain six comparisons in seven scenarios")
assert(nrow(scenario_estimands) == 7L && nrow(decision) == 1L,
       "Confirmatory estimand outputs have the wrong dimensions")
assert(nrow(sensitivity) == 4L && sum(sensitivity$role == "primary") == 1L,
       "Prediction-margin sensitivity output is incomplete")

## Independently recompute the locked primary estimands and bounds from raw
## replication-level rows, rather than trusting the analysis output.
scenario_support <- scenario_support_se <- numeric(nrow(SCENARIOS))
scenario_prediction <- scenario_prediction_se <- numeric(nrow(SCENARIOS))
paired_counts <- integer(nrow(SCENARIOS))
for (i in seq_len(nrow(SCENARIOS))) {
  dat <- all_primary[[SCENARIOS$scenario_id[i]]]
  lhs <- dat[dat$method == "PanIC-CF" & dat$method_failed == 0L, ]
  rhs <- dat[
    dat$method == CONFIG$primary_cv_method & dat$method_failed == 0L,
  ]
  pair <- merge(lhs, rhs, by = "replication", suffixes = c("_lhs", "_rhs"))
  paired_counts[i] <- nrow(pair)
  support_delta <- pair$wrong_lhs - pair$wrong_rhs
  prediction_delta <-
    (pair$test_deviance_lhs - pair$test_deviance_rhs) /
    pair$test_deviance_rhs
  scenario_support[i] <- mean(support_delta)
  scenario_support_se[i] <- mcse(support_delta)
  scenario_prediction[i] <- mean(prediction_delta)
  scenario_prediction_se[i] <- mcse(prediction_delta)
}
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
       "Support decision was not reproduced from raw rows")
assert(close(decision$prediction_estimate, prediction_estimate) &&
         close(decision$prediction_mcse, prediction_mcse) &&
         close(decision$prediction_upper_one_sided_95, prediction_upper),
       "Prediction decision was not reproduced from raw rows")
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
         ), "The locked joint decision rule was applied incorrectly")

required_artifacts <- c(
  "diagnostic_summary.csv", "calibration_target_summary.csv",
  "manuscript_generated/table_primary_support.tex",
  "manuscript_generated/table_confirmatory_decision.tex",
  "sessionInfo.txt", "environment.txt", "locked_configuration.rds"
)
assert(all(file.exists(file.path(results_dir, required_artifacts))),
       "One or more required confirmatory artifacts are missing")
support_table_text <- paste(readLines(
  file.path(results_dir, "manuscript_generated", "table_primary_support.tex"),
  warn = FALSE
), collapse = "\n")
assert(!grepl("BIC-active \\(exploratory\\)", support_table_text) &&
         grepl("BIC-like", support_table_text),
       "Reader-facing BIC-like labels were not normalized")
assert(grepl(
  "exploratory active-count analogues outside the scope of the Gaussian",
  support_table_text, fixed = TRUE
), "The non-Gaussian BIC-like qualification is missing")
assert(!grepl(
  "second-confirmation|revised PanIC-CF", support_table_text,
  ignore.case = TRUE
), "A manuscript table exposes internal process wording")

cat("PASS: confirmatory results and decision rule are verified.\n")
cat("Joint claim decision:", decision$joint_claim_pass, "\n")
cat("Support upper bound:", format(support_upper, digits = 8), "\n")
cat("Prediction upper bound:", format(prediction_upper, digits = 8), "\n")
