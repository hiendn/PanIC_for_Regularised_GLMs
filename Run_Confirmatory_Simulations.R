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
  normalizePath("Run_Confirmatory_Simulations.R")
}
candidate_dir <- dirname(script_path)
.libPaths(c(file.path(candidate_dir, "Rlib"), .libPaths()))
source(file.path(candidate_dir, "Method_Lock_Verification.R"))
verify_method_lock(candidate_dir)
source(file.path(candidate_dir, "Simulation_Config.R"))
PRODUCTION_MASTER_SEED <- CONFIG$master_seed
smoke_mode <- "--smoke" %in% args
if (smoke_mode) CONFIG$master_seed <- CONFIG$smoke_master_seed
source(file.path(candidate_dir, "PanIC_CF_Functions.R"))

n_rep <- as.integer(arg_value(
  "--n-rep", if (smoke_mode) 2L else CONFIG$n_rep
))
cores <- max(1L, as.integer(arg_value("--cores", CONFIG$requested_cores)))
scenario_filter <- arg_value("--scenarios", "all")
output_name <- arg_value("--output-dir", CONFIG$default_output_dir)
overwrite <- "--overwrite" %in% args
if (!is.finite(n_rep) || n_rep < 1L) {
  stop("--n-rep must be a positive integer")
}
results_dir <- if (grepl("^/", output_name)) {
  output_name
} else {
  file.path(candidate_dir, output_name)
}
is_locked_output <- identical(
  normalizePath(dirname(results_dir), mustWork = FALSE),
  normalizePath(candidate_dir, mustWork = TRUE)
) && identical(basename(results_dir), CONFIG$default_output_dir)
if (smoke_mode && is_locked_output) {
  stop("Smoke mode must use an output directory other than results")
}
if (!smoke_mode && n_rep != CONFIG$n_rep) {
  stop(
    "The production seed requires exactly ", CONFIG$n_rep,
    " replications; add --smoke and use a different output directory for ",
    "any diagnostic run"
  )
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

if (scenario_filter == "all") {
  scenario_indices <- seq_len(nrow(SCENARIOS))
} else {
  requested <- strsplit(scenario_filter, ",", fixed = TRUE)[[1L]]
  scenario_indices <- match(requested, SCENARIOS$scenario_id)
  if (anyNA(scenario_indices)) stop("Unknown scenario id in --scenarios")
}
if (!smoke_mode && scenario_filter != "all") {
  stop("The locked production run must execute all seven scenarios together")
}

configuration <- data.frame(
  item = c(
    "method_version", "seed_role", "master_seed",
    "production_master_seed", "replications_requested",
    "independent_test_size", "dimension", "active_slopes", "radius_min",
    "radius_max", "radius_points", "kappa_grid",
    "calibration_default_kappa", "calibration_half_splits",
    "calibration_directions_per_split", "calibration_rows",
    "calibration_weight", "cv_folds", "cv_primary_rule",
    "support_definition", "primary_support_estimand",
    "primary_support_decision", "prediction_guardrail_estimand",
    "prediction_noninferiority_margin",
    "prediction_noninferiority_decision", "one_sided_confidence_level",
    "prediction_sensitivity_margins", "workers_requested",
    "glmnet_nlambda", "glmnet_lambda_min_ratio", "glmnet_threshold",
    "glmnet_maxit", "interpolation_radius_tolerance",
    "pilot_coefficient_cap"
  ),
  value = c(
    "PanIC-CF confirmatory simulation",
    if (smoke_mode) "dedicated-smoke" else "production",
    CONFIG$master_seed, PRODUCTION_MASTER_SEED, n_rep,
    CONFIG$n_test, CONFIG$d, length(CONFIG$active), CONFIG$radius_min,
    CONFIG$radius_max, CONFIG$radius_points,
    "31 log-spaced values from 0.01 to 100",
    CONFIG$calibration_default_kappa, CONFIG$calibration_half_splits,
    CONFIG$calibration_directions,
    CONFIG$calibration_half_splits * CONFIG$calibration_directions,
    "sqrt(log(log(n_validation + exp(exp(1)))))", CONFIG$cv_folds,
    "minimum mean five-fold validation loss on radius grid",
    "beta_hat != 0 (literal fitted nonzero)",
    "equal-weight seven-scenario mean paired (FP+FN)_PanIC-(FP+FN)_CV",
    "one-sided 95% upper Monte Carlo bound < 0",
    paste0(
      "equal-weight seven-scenario mean paired relative deviance change ",
      "(D_PanIC-D_CV)/D_CV"
    ),
    CONFIG$prediction_noninferiority_margin,
    paste0(
      "one-sided 95% upper Monte Carlo bound < ",
      CONFIG$prediction_noninferiority_margin
    ),
    1 - CONFIG$primary_support_alpha,
    paste(CONFIG$prediction_sensitivity_margins, collapse = ";"), cores,
    CONFIG$glmnet_nlambda, CONFIG$glmnet_lambda_min_ratio,
    CONFIG$glmnet_threshold, CONFIG$glmnet_maxit,
    CONFIG$interpolation_radius_tolerance, CONFIG$pilot_coefficient_cap
  ),
  stringsAsFactors = FALSE
)

existing_lock <- file.path(results_dir, "locked_configuration.rds")
if (file.exists(existing_lock) && !overwrite) {
  prior_lock <- readRDS(existing_lock)
  if (!isTRUE(all.equal(prior_lock, CONFIG, check.attributes = TRUE))) {
    stop("Existing output directory has a different locked configuration")
  }
}
write.csv(configuration, file.path(results_dir, "configuration.csv"),
          row.names = FALSE)
write.csv(SCENARIOS, file.path(results_dir, "scenario_registry.csv"),
          row.names = FALSE)
saveRDS(CONFIG, existing_lock, compress = "xz")

run_scenario <- function(i) {
  scenario <- SCENARIOS[i, , drop = FALSE]
  id <- scenario$scenario_id[[1L]]
  paths <- file.path(
    results_dir,
    paste0(id, c("_primary.csv", "_diagnostics.csv",
                 "_calibration_rows.csv"))
  )
  raw_path <- file.path(results_dir, paste0(id, "_raw.rds"))
  if (!overwrite && all(file.exists(c(paths, raw_path)))) {
    message("Skipping completed scenario: ", id)
    return(invisible(NULL))
  }
  message(
    "Running ", id, " with ", n_rep, " attempted replications on ",
    cores, " workers"
  )
  worker <- function(r) {
    tryCatch(
      run_candidate_replication(scenario, i, r),
      error = function(e) unexpected_candidate_error_result(
        scenario, i, r, conditionMessage(e), CONFIG
      )
    )
  }
  started <- proc.time()[["elapsed"]]
  if (.Platform$OS.type == "unix" && cores > 1L) {
    pieces <- parallel::mclapply(
      seq_len(n_rep), worker, mc.cores = cores, mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    pieces <- lapply(seq_len(n_rep), worker)
  }
  elapsed <- proc.time()[["elapsed"]] - started
  primary <- do.call(rbind, lapply(pieces, `[[`, "primary"))
  diagnostic <- do.call(rbind, lapply(pieces, `[[`, "diagnostic"))
  calibration <- do.call(rbind, lapply(pieces, `[[`, "calibration"))
  write.csv(primary, paths[[1L]], row.names = FALSE)
  write.csv(diagnostic, paths[[2L]], row.names = FALSE)
  write.csv(calibration, paths[[3L]], row.names = FALSE)
  saveRDS(pieces, raw_path, compress = "xz")
  writeLines(
    c(
      paste0("scenario_id=", id),
      paste0("attempted_replications=", n_rep),
      paste0("workers=", cores),
      paste0("wall_elapsed_seconds=", sprintf("%.6f", elapsed)),
      paste0("full_path_failures=", sum(diagnostic$full_path_failed)),
      paste0("calibration_failures=", sum(diagnostic$calibration_failed)),
      paste0("cv_failures=", sum(diagnostic$cv_failed))
    ),
    file.path(results_dir, paste0(id, "_runtime.txt"))
  )
  message("Completed ", id, " in ", sprintf("%.1f", elapsed), " seconds")
}

started_all <- proc.time()[["elapsed"]]
for (i in scenario_indices) run_scenario(i)
elapsed_all <- proc.time()[["elapsed"]] - started_all

diagnostic_files <- file.path(
  results_dir, paste0(SCENARIOS$scenario_id, "_diagnostics.csv")
)
diagnostic_files <- diagnostic_files[file.exists(diagnostic_files)]
if (length(diagnostic_files)) {
  diagnostics <- do.call(rbind, lapply(
    diagnostic_files, read.csv, stringsAsFactors = FALSE
  ))
  split_seed_columns <- paste0(
    "calibration_split_seed_", seq_len(CONFIG$calibration_half_splits)
  )
  split_seed_matrix <- t(vapply(
    seq_len(nrow(diagnostics)),
    function(j) {
      scenario_index <- match(
        diagnostics$scenario_id[j], SCENARIOS$scenario_id
      )
      seed_streams(
        scenario_index, diagnostics$replication[j]
      )$calibration_splits
    },
    integer(CONFIG$calibration_half_splits)
  ))
  colnames(split_seed_matrix) <- split_seed_columns
  seed_ledger <- cbind(
    diagnostics[c(
      "scenario_id", "replication", "training_seed", "cv_fold_seed",
      "test_seed"
    )],
    as.data.frame(split_seed_matrix)
  )
  seed_ledger <- seed_ledger[order(
    match(seed_ledger$scenario_id, SCENARIOS$scenario_id),
    seed_ledger$replication
  ), ]
  write.csv(seed_ledger, file.path(results_dir, "seed_ledger.csv"),
            row.names = FALSE)
}

writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo.txt"))
method_lock_hash <- system2(
  "shasum", c(
    "-a", "256",
    shQuote(file.path(candidate_dir, "METHOD_LOCK_SHA256.txt"))
  ),
  stdout = TRUE
)
writeLines(
  c(
    paste0("R=", R.version.string),
    paste0("platform=", R.version$platform),
    paste0("glmnet=", as.character(packageVersion("glmnet"))),
    paste0("seed_role=", if (smoke_mode) "dedicated-smoke" else "production"),
    paste0("master_seed=", CONFIG$master_seed),
    paste0("production_master_seed=", PRODUCTION_MASTER_SEED),
    paste0("method_lock_sha256=", substr(method_lock_hash[[1L]], 1L, 64L)),
    paste0("overall_wall_elapsed_seconds=", sprintf("%.6f", elapsed_all)),
    paste0("completed_at=", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ),
  file.path(results_dir, "environment.txt")
)
message("All requested scenarios completed in ",
        sprintf("%.1f", elapsed_all), " seconds")
