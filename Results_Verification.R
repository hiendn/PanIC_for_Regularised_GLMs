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
} else normalizePath("Results_Verification.R")
production_dir <- dirname(script_path)
source(file.path(production_dir, "Simulation_Config.R"))
output_name <- arg_value("--output-dir", "results")
results_dir <- if (grepl("^/", output_name)) output_name else
  file.path(production_dir, output_name)

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}
configuration <- read.csv(file.path(results_dir, "configuration.csv"),
                          stringsAsFactors = FALSE)
n_rep <- as.integer(configuration$value[
  configuration$item == "replications_requested"
])
assert(length(n_rep) == 1L && is.finite(n_rep),
       "Cannot recover the requested replication count")
validation <- read.csv(file.path(results_dir, "verification_checks.csv"),
                       stringsAsFactors = FALSE)
assert(all(validation$passed == 1L),
       "A deterministic implementation check failed")

all_diagnostics <- list()
for (id in SCENARIOS$scenario_id) {
  primary_path <- file.path(results_dir, paste0(id, "_primary.csv"))
  diagnostic_path <- file.path(results_dir, paste0(id, "_diagnostics.csv"))
  calibration_path <- file.path(results_dir, paste0(id, "_calibration_rows.csv"))
  assert(file.exists(primary_path), paste0("Missing primary file for ", id))
  assert(file.exists(diagnostic_path), paste0("Missing diagnostics file for ", id))
  assert(file.exists(calibration_path), paste0("Missing calibration file for ", id))
  primary <- read.csv(primary_path, stringsAsFactors = FALSE)
  diagnostic <- read.csv(diagnostic_path, stringsAsFactors = FALSE)
  calibration <- read.csv(calibration_path, stringsAsFactors = FALSE)
  all_diagnostics[[id]] <- diagnostic
  assert(nrow(primary) == 3L * n_rep,
         paste0(id, ": expected three methods per attempted replication"))
  assert(nrow(diagnostic) == n_rep,
         paste0(id, ": diagnostics do not cover all attempted replications"))
  assert(nrow(calibration) == 10L * n_rep,
         paste0(id, ": expected ten calibration rows per replication"))
  assert(identical(sort(unique(diagnostic$replication)), seq_len(n_rep)),
         paste0(id, ": incomplete replication indices"))
  assert(all(table(primary$replication) == 3L),
         paste0(id, ": method row count differs by replication"))
  assert(all(table(calibration$replication) == 10L),
         paste0(id, ": calibration row count differs by replication"))
  assert(all(table(calibration$replication, calibration$split_repeat) == 2L),
         paste0(id, ": a half-split is not scored in both directions"))
  assert(all(table(calibration$replication, calibration$direction) == 5L),
         paste0(id, ": a direction is not represented five times"))
  ok_target <- is.finite(calibration$raw_cross_signed_radius)
  if (any(ok_target)) {
    expected_projection <- pmin(
      CONFIG$radius_max,
      pmax(CONFIG$radius_min,
           calibration$raw_cross_signed_radius[ok_target])
    )
    assert(max(abs(
      expected_projection -
        calibration$projected_cross_signed_radius[ok_target]
    )) < 1e-12, paste0(id, ": target projection mismatch"))
    expected_weight <- log(log(
      calibration$n_validation[ok_target] + exp(exp(1))
    ))
    assert(max(abs(expected_weight - calibration$weight[ok_target])) < 1e-12,
           paste0(id, ": calibration weight mismatch"))
  }
  panic <- primary[primary$method == "PanIC-CF", ]
  cv <- primary[primary$method == "CV", ]
  assert(all(panic$method_failed == diagnostic$full_path_failed),
         paste0(id, ": PanIC-CF numerical-failure flag is inconsistent"))
  assert(all(panic$calibration_failed == diagnostic$calibration_failed),
         paste0(id, ": PanIC-CF calibration-failure flag is inconsistent"))
  assert(all(panic$default_used == diagnostic$default_used),
         paste0(id, ": PanIC-CF default-use flag is inconsistent"))
  assert(all(
    diagnostic$default_used ==
      as.integer(diagnostic$calibration_failed == 1L &
                   diagnostic$full_path_failed == 0L)
  ), paste0(id, ": prespecified default-use rule is inconsistent"))
  assert(all(cv$method_failed == diagnostic$cv_failed),
         paste0(id, ": CV failure flag is inconsistent"))
  good_kappa <- diagnostic$full_path_failed == 0L
  kappa_on_grid <- vapply(
    diagnostic$kappa_hat[good_kappa],
    function(value) min(abs(value - CONFIG$kappa_grid)) < 1e-12,
    logical(1)
  )
  assert(all(kappa_on_grid),
         paste0(id, ": selected kappa is off-grid"))
  finite_radius_error <- diagnostic$maximum_active_radius_error[
    is.finite(diagnostic$maximum_active_radius_error)
  ]
  assert(!length(finite_radius_error) ||
           max(finite_radius_error) <= CONFIG$interpolation_radius_tolerance,
         paste0(id, ": interpolation tolerance exceeded"))
}

diagnostics <- do.call(rbind, all_diagnostics)
seed_ledger <- read.csv(file.path(results_dir, "seed_ledger.csv"),
                        stringsAsFactors = FALSE)
assert(nrow(seed_ledger) == nrow(SCENARIOS) * n_rep,
       "Seed ledger does not have one row per attempted replication")
assert(!anyDuplicated(seed_ledger[c("scenario_id", "replication")]),
       "Seed ledger duplicates a scenario/replication pair")
seed_columns <- setdiff(names(seed_ledger), c("scenario_id", "replication"))
assert(all(apply(seed_ledger[seed_columns], 1L, function(x) !anyDuplicated(x))),
       "Seed streams collide within a replication")
assert(all(diagnostics$calibration_rows_expected == 10L),
       "A diagnostic row reports a non-locked calibration-row count")

required_summaries <- c(
  "simulation_summary.csv", "paired_panic_cf_vs_cv.csv",
  "diagnostic_summary.csv", "calibration_target_summary.csv",
  "table_primary_support.tex",
  "table_primary_performance.tex",
  "table_calibration.tex",
  "sessionInfo.txt", "environment.txt"
)
assert(all(file.exists(file.path(results_dir, required_summaries))),
       "One or more summary artifacts are missing")

cat("All confirmatory-result integrity checks passed.\n")
cat("Recorded full-path failures:", sum(diagnostics$full_path_failed), "\n")
cat("Recorded calibration failures:", sum(diagnostics$calibration_failed), "\n")
cat("Recorded CV failures:", sum(diagnostics$cv_failed), "\n")
