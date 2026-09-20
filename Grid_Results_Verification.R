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
  "PanIC-CF", "BIC-like", CONFIG$primary_cv_method
)
calibration_by_grid <- list()
primary_by_grid <- list()
diagnostic_by_grid <- list()
for (m in CONFIG$grid_sensitivity_points) {
  primary <- read_result(paste0("grid_sensitivity_m", m, "_primary.csv"))
  diagnostic <- read_result(
    paste0("grid_sensitivity_m", m, "_diagnostics.csv")
  )
  calibration <- read_result(
    paste0("grid_sensitivity_m", m, "_calibration_rows.csv")
  )
  calibration_by_grid[[as.character(m)]] <- calibration
  primary_by_grid[[as.character(m)]] <- primary
  diagnostic_by_grid[[as.character(m)]] <- diagnostic
  assert(nrow(primary) == 3L * n_rep,
         paste0("m=", m, ": primary row count mismatch"))
  assert(nrow(diagnostic) == n_rep,
         paste0("m=", m, ": diagnostic row count mismatch"))
  assert(nrow(calibration) == 10L * n_rep,
         paste0("m=", m, ": calibration row count mismatch"))
  assert(all(vapply(
    split(primary$method, primary$replication),
    function(methods) setequal(methods, method_order), logical(1)
  )), paste0("m=", m, ": method set mismatch"))
  expected_weight <- sqrt(log(log(
    calibration$n_validation + exp(exp(1))
  )))
  target_ok <- is.finite(calibration$raw_cross_signed_radius)
  assert(max(abs(
    calibration$weight[target_ok] - expected_weight[target_ok]
  )) < 1e-12, paste0("m=", m, ": calibration weight mismatch"))
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
    "weight"
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
assert(nrow(summary) == 9L,
       "Grid summary must contain three methods on three grids")
assert(nrow(grid_contrasts) == 9L,
       "Grid paired table must contain three methods and three grid pairs")
assert(nrow(method_contrasts) == 9L,
       "Grid method table must contain three comparisons on three grids")
assert(all(summary$attempted_replications == n_rep),
       "Grid summary replication count mismatch")

## Independently reconstruct every summary field displayed in the grid table
## from the replication-level primary and diagnostic rows.
all_grid_primary <- do.call(rbind, primary_by_grid)
all_grid_diagnostic <- do.call(rbind, diagnostic_by_grid)
grid_mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x) / sqrt(length(x))
}
grid_metrics <- c(
  "fp", "fn", "fpr", "fnr", "exact", "wrong",
  "signed_attained_radius_error", "test_deviance",
  "selected_grid_radius", "attained_radius",
  "selected_lower_endpoint", "selected_upper_endpoint"
)
expected_grid_rows <- lapply(
  split(
    all_grid_primary,
    interaction(
      all_grid_primary$radius_points, all_grid_primary$method, drop = TRUE
    )
  ),
  function(dat) {
    out <- dat[1L, c("radius_points", "method")]
    out$grid_spacing <- (CONFIG$radius_max - CONFIG$radius_min) /
      (out$radius_points - 1L)
    out$attempted_replications <- nrow(dat)
    out$method_failures <- sum(dat$method_failed)
    for (metric in grid_metrics) {
      out[[paste0(metric, "_mean")]] <- mean(dat[[metric]], na.rm = TRUE)
      out[[paste0(metric, "_mcse")]] <- grid_mcse(dat[[metric]])
    }
    diagnostic <- all_grid_diagnostic[
      all_grid_diagnostic$radius_points == out$radius_points,
    ]
    out$full_path_failures <- sum(diagnostic$full_path_failed)
    out$calibration_failures <- sum(diagnostic$calibration_failed)
    out$calibration_default_uses <- sum(diagnostic$default_used)
    out$cv_failures <- sum(diagnostic$cv_failed)
    out$mean_shared_runtime <- mean(
      diagnostic$elapsed_seconds, na.rm = TRUE
    )
    out$mcse_shared_runtime <- grid_mcse(diagnostic$elapsed_seconds)
    out$maximum_radius_interpolation_error <- max(
      diagnostic$maximum_active_radius_error, na.rm = TRUE
    )
    out
  }
)
expected_grid_summary <- do.call(rbind, expected_grid_rows)
reported_grid_key <- paste(summary$radius_points, summary$method, sep = "\r")
expected_grid_key <- paste(
  expected_grid_summary$radius_points,
  expected_grid_summary$method,
  sep = "\r"
)
grid_match <- match(reported_grid_key, expected_grid_key)
assert(
  !anyNA(grid_match) && !anyDuplicated(reported_grid_key) &&
    !anyDuplicated(expected_grid_key),
  "Grid-summary radius/method keys are incomplete or duplicated"
)
expected_grid_summary <- expected_grid_summary[grid_match, , drop = FALSE]
rownames(expected_grid_summary) <- NULL
rownames(summary) <- NULL
grid_character_columns <- "method"
grid_numeric_columns <- setdiff(names(summary), grid_character_columns)
assert(
  identical(
    summary[grid_character_columns],
    expected_grid_summary[grid_character_columns]
  ) &&
    isTRUE(all.equal(
      summary[grid_numeric_columns],
      expected_grid_summary[grid_numeric_columns],
      tolerance = 1e-12, check.attributes = FALSE
    )),
  "Grid table summary was not reproduced from raw rows"
)
grid_table_path <- file.path(results_dir, "table_grid_sensitivity.tex")
assert(file.exists(grid_table_path),
       "Missing canonical grid-sensitivity table mirror")
assert(file.exists(file.path(results_dir, "grid_locked_configuration.rds")),
       "Missing locked grid configuration")
expected_grid_table_bytes <- panic_table_bytes(
  render_grid_sensitivity_table(summary)
)
observed_grid_table_size <- file.info(grid_table_path)$size
observed_grid_table_bytes <- readBin(
  grid_table_path, what = "raw", n = observed_grid_table_size
)
assert(
  identical(
    observed_grid_table_size,
    as.numeric(length(expected_grid_table_bytes))
  ) && identical(observed_grid_table_bytes, expected_grid_table_bytes),
  "Canonical grid-sensitivity table mirror is stale or non-deterministic"
)
grid_table_text <- paste(readLines(
  grid_table_path,
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
