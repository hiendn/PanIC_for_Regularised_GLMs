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
  normalizePath("Grid_Sensitivity_Study.R")
}
candidate_dir <- dirname(script_path)
.libPaths(c(file.path(candidate_dir, "Rlib"), .libPaths()))
source(file.path(candidate_dir, "Method_Lock_Verification.R"))
verify_method_lock(candidate_dir)
source(file.path(candidate_dir, "Simulation_Config.R"))
source(file.path(candidate_dir, "Manuscript_Table_Rendering.R"))

PRODUCTION_CONFIRMATORY_MASTER_SEED <- CONFIG$master_seed
PRODUCTION_SENSITIVITY_MASTER_SEED <- CONFIG$grid_sensitivity_master_seed
smoke_mode <- "--smoke" %in% args
CONFIRMATORY_MASTER_SEED <- if (smoke_mode) {
  CONFIG$smoke_master_seed
} else {
  PRODUCTION_CONFIRMATORY_MASTER_SEED
}
SENSITIVITY_MASTER_SEED <- if (smoke_mode) {
  CONFIG$smoke_grid_sensitivity_master_seed
} else {
  PRODUCTION_SENSITIVITY_MASTER_SEED
}
CONFIG$master_seed <- SENSITIVITY_MASTER_SEED
source(file.path(candidate_dir, "PanIC_CF_Functions.R"))

n_rep <- as.integer(arg_value(
  "--n-rep", if (smoke_mode) 2L else CONFIG$grid_sensitivity_replications
))
cores <- max(1L, as.integer(arg_value("--cores", CONFIG$requested_cores)))
output_name <- arg_value("--output-dir", CONFIG$default_output_dir)
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
if (!smoke_mode && n_rep != CONFIG$grid_sensitivity_replications) {
  stop(
    "The production grid seed requires exactly ",
    CONFIG$grid_sensitivity_replications,
    " replications; add --smoke and use a different output directory for ",
    "any diagnostic run"
  )
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

grid_sizes <- CONFIG$grid_sensitivity_points
scenario <- data.frame(
  scenario_id = "linear_iid_n1000_grid_sensitivity",
  family = "gaussian", n = 1000L, rho = 0,
  stringsAsFactors = FALSE
)
scenario_seed_index <- 1L

for (m in grid_sizes) {
  message("Running second-confirmation grid sensitivity m=", m,
          " with common random numbers")
  worker <- function(replication) {
    tryCatch(
      run_candidate_replication(
        scenario, scenario_seed_index, replication,
        radii = radius_grid(m), config = CONFIG
      ),
      error = function(e) unexpected_candidate_error_result(
        scenario, scenario_seed_index, replication,
        conditionMessage(e), CONFIG
      )
    )
  }
  started <- proc.time()[["elapsed"]]
  pieces <- if (.Platform$OS.type == "unix" && cores > 1L) {
    parallel::mclapply(
      seq_len(n_rep), worker, mc.cores = cores, mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(seq_len(n_rep), worker)
  }
  wall_elapsed <- proc.time()[["elapsed"]] - started
  primary <- do.call(rbind, lapply(pieces, `[[`, "primary"))
  diagnostic <- do.call(rbind, lapply(pieces, `[[`, "diagnostic"))
  calibration <- do.call(rbind, lapply(pieces, `[[`, "calibration"))
  primary$radius_points <- m
  diagnostic$radius_points <- m
  calibration$radius_points <- m
  stem <- file.path(results_dir, paste0("grid_sensitivity_m", m))
  write.csv(primary, paste0(stem, "_primary.csv"), row.names = FALSE)
  write.csv(diagnostic, paste0(stem, "_diagnostics.csv"), row.names = FALSE)
  write.csv(calibration, paste0(stem, "_calibration_rows.csv"),
            row.names = FALSE)
  saveRDS(pieces, paste0(stem, "_raw.rds"), compress = "xz")
  writeLines(
    c(
      paste0("radius_points=", m),
      paste0("attempted_replications=", n_rep),
      paste0("workers=", cores),
      paste0("wall_elapsed_seconds=", sprintf("%.6f", wall_elapsed)),
      paste0("sensitivity_master_seed=", SENSITIVITY_MASTER_SEED),
      paste0("confirmatory_master_seed=", CONFIRMATORY_MASTER_SEED),
      paste0("full_path_failures=", sum(diagnostic$full_path_failed)),
      paste0("calibration_failures=", sum(diagnostic$calibration_failed)),
      paste0("calibration_default_uses=", sum(diagnostic$default_used)),
      paste0("cv_failures=", sum(diagnostic$cv_failed))
    ),
    paste0(stem, "_runtime.txt")
  )
  message("Completed grid m=", m, " in ",
          sprintf("%.1f", wall_elapsed), " seconds")
}

read_grid <- function(suffix) {
  do.call(rbind, lapply(grid_sizes, function(m) {
    read.csv(
      file.path(results_dir, paste0("grid_sensitivity_m", m, suffix)),
      stringsAsFactors = FALSE
    )
  }))
}
primary <- read_grid("_primary.csv")
diagnostic <- read_grid("_diagnostics.csv")
calibration <- read_grid("_calibration_rows.csv")
mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

metrics <- c(
  "fp", "fn", "fpr", "fnr", "exact", "wrong",
  "signed_attained_radius_error", "test_deviance",
  "selected_grid_radius", "attained_radius",
  "selected_lower_endpoint", "selected_upper_endpoint"
)
summary_rows <- lapply(
  split(primary, interaction(primary$radius_points, primary$method,
                             drop = TRUE)),
  function(dat) {
    out <- dat[1L, c("radius_points", "method")]
    out$grid_spacing <- (CONFIG$radius_max - CONFIG$radius_min) /
      (out$radius_points - 1L)
    out$attempted_replications <- nrow(dat)
    out$method_failures <- sum(dat$method_failed)
    for (metric in metrics) {
      out[[paste0(metric, "_mean")]] <- mean(dat[[metric]], na.rm = TRUE)
      out[[paste0(metric, "_mcse")]] <- mcse(dat[[metric]])
    }
    d <- diagnostic[diagnostic$radius_points == out$radius_points, ]
    out$full_path_failures <- sum(d$full_path_failed)
    out$calibration_failures <- sum(d$calibration_failed)
    out$calibration_default_uses <- sum(d$default_used)
    out$cv_failures <- sum(d$cv_failed)
    out$mean_shared_runtime <- mean(d$elapsed_seconds, na.rm = TRUE)
    out$mcse_shared_runtime <- mcse(d$elapsed_seconds)
    out$maximum_radius_interpolation_error <- max(
      d$maximum_active_radius_error, na.rm = TRUE
    )
    out
  }
)
grid_summary <- do.call(rbind, summary_rows)
method_order <- c(
  "PanIC-CF", CONFIG$original_sensitivity_method, "BIC-like",
  CONFIG$primary_cv_method, CONFIG$secondary_cv_method
)
grid_summary <- grid_summary[order(
  grid_summary$radius_points, match(grid_summary$method, method_order)
), ]
write.csv(grid_summary,
          file.path(results_dir, "grid_sensitivity_summary.csv"),
          row.names = FALSE)

## Paired differences across grid sizes exploit common random numbers.
paired_rows <- list()
cursor <- 1L
grid_pairs <- list(c(61L, 121L), c(121L, 241L), c(61L, 241L))
contrast_metrics <- c(
  "wrong", "exact", "signed_attained_radius_error", "test_deviance",
  "selected_grid_radius", "attained_radius"
)
for (method in method_order) {
  method_data <- primary[
    primary$method == method & primary$method_failed == 0L,
  ]
  for (pair in grid_pairs) {
    lower <- method_data[method_data$radius_points == pair[1L], ]
    upper <- method_data[method_data$radius_points == pair[2L], ]
    paired <- merge(lower, upper, by = "replication",
                    suffixes = c("_lower", "_upper"))
    out <- data.frame(
      method = method, lower_radius_points = pair[1L],
      upper_radius_points = pair[2L], paired_replications = nrow(paired),
      stringsAsFactors = FALSE
    )
    for (metric in contrast_metrics) {
      delta <- paired[[paste0(metric, "_upper")]] -
        paired[[paste0(metric, "_lower")]]
      out[[paste0(metric, "_difference")]] <- mean(delta)
      out[[paste0(metric, "_difference_mcse")]] <- mcse(delta)
    }
    paired_rows[[cursor]] <- out
    cursor <- cursor + 1L
  }
}
paired_grid_contrasts <- do.call(rbind, paired_rows)
write.csv(
  paired_grid_contrasts,
  file.path(results_dir, "grid_sensitivity_paired_contrasts.csv"),
  row.names = FALSE
)

## Within-grid method comparisons preserve the revised-vs-CV and
## revised-vs-original distinctions used by the main study.
method_pairs <- list(
  c("PanIC-CF", CONFIG$primary_cv_method),
  c(CONFIG$original_sensitivity_method, CONFIG$primary_cv_method),
  c("PanIC-CF", CONFIG$original_sensitivity_method),
  c("PanIC-CF", CONFIG$secondary_cv_method)
)
method_rows <- list()
cursor <- 1L
for (m in grid_sizes) {
  grid_data <- primary[primary$radius_points == m, ]
  for (pair in method_pairs) {
    lhs <- grid_data[
      grid_data$method == pair[1L] & grid_data$method_failed == 0L,
    ]
    rhs <- grid_data[
      grid_data$method == pair[2L] & grid_data$method_failed == 0L,
    ]
    paired <- merge(lhs, rhs, by = "replication",
                    suffixes = c("_lhs", "_rhs"))
    out <- data.frame(
      radius_points = m, lhs_method = pair[1L], rhs_method = pair[2L],
      comparison = paste(pair[1L], "minus", pair[2L]),
      paired_replications = nrow(paired), stringsAsFactors = FALSE
    )
    for (metric in contrast_metrics) {
      delta <- paired[[paste0(metric, "_lhs")]] -
        paired[[paste0(metric, "_rhs")]]
      out[[paste0(metric, "_difference")]] <- mean(delta)
      out[[paste0(metric, "_difference_mcse")]] <- mcse(delta)
    }
    method_rows[[cursor]] <- out
    cursor <- cursor + 1L
  }
}
method_contrasts <- do.call(rbind, method_rows)
write.csv(
  method_contrasts,
  file.path(results_dir, "grid_sensitivity_method_contrasts.csv"),
  row.names = FALSE
)

## Verify that data, seeds, split targets, and both weights are identical
## across the three grids before any grid comparison is retained.
calibration_key <- c("replication", "split_repeat", "direction")
common_columns <- c(
  calibration_key, "training_seed", "split_seed", "test_seed",
  "raw_cross_signed_radius", "weight", "original_weight"
)
reference <- calibration[
  calibration$radius_points == 121L, common_columns
]
for (m in c(61L, 241L)) {
  candidate <- calibration[
    calibration$radius_points == m, common_columns
  ]
  comparison <- merge(
    reference, candidate, by = calibration_key,
    suffixes = c("_reference", "_candidate")
  )
  exact_columns <- c("training_seed", "split_seed", "test_seed")
  exact_ok <- all(vapply(exact_columns, function(name) {
    all(comparison[[paste0(name, "_reference")]] ==
          comparison[[paste0(name, "_candidate")]])
  }, logical(1)))
  numeric_columns <- c(
    "raw_cross_signed_radius", "weight", "original_weight"
  )
  numeric_ok <- all(vapply(numeric_columns, function(name) {
    max(abs(
      comparison[[paste0(name, "_reference")]] -
        comparison[[paste0(name, "_candidate")]]
    ), na.rm = TRUE) <= 1e-12
  }, logical(1)))
  if (nrow(comparison) != 10L * n_rep || !exact_ok || !numeric_ok) {
    stop("Common-random-number grid comparison failed for m=", m)
  }
}

grid_configuration <- data.frame(
  item = c(
    "seed_role", "sensitivity_master_seed", "confirmatory_master_seed",
    "production_sensitivity_master_seed",
    "production_confirmatory_master_seed",
    "scenario_seed_index", "replications", "workers", "radius_points",
    "radius_interval", "methods"
  ),
  value = c(
    if (smoke_mode) "dedicated-smoke" else "production",
    SENSITIVITY_MASTER_SEED, CONFIRMATORY_MASTER_SEED,
    PRODUCTION_SENSITIVITY_MASTER_SEED,
    PRODUCTION_CONFIRMATORY_MASTER_SEED,
    scenario_seed_index, n_rep, cores, paste(grid_sizes, collapse = ";"),
    paste0("[", CONFIG$radius_min, ",", CONFIG$radius_max, "]"),
    paste(method_order, collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write.csv(
  grid_configuration,
  file.path(results_dir, "grid_sensitivity_configuration.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    seed_role = if (smoke_mode) "dedicated-smoke" else "production",
    confirmatory_master_seed = CONFIRMATORY_MASTER_SEED,
    sensitivity_master_seed = SENSITIVITY_MASTER_SEED,
    production_confirmatory_master_seed =
      PRODUCTION_CONFIRMATORY_MASTER_SEED,
    production_sensitivity_master_seed =
      PRODUCTION_SENSITIVITY_MASTER_SEED,
    replications = n_rep, grid_sizes = grid_sizes,
    scenario = scenario, config = CONFIG
  ),
  file.path(results_dir, "grid_locked_configuration.rds"), compress = "xz"
)
## Re-read the retained summary before presentation rendering so the CSV is
## the single deterministic formatting contract at decimal-rounding ties.
render_grid_summary <- read.csv(
  file.path(results_dir, "grid_sensitivity_summary.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
grid_tables <- setNames(
  list(render_grid_sensitivity_table(render_grid_summary)),
  "table_grid_sensitivity.tex"
)
write_manuscript_table_subset(grid_tables, results_dir)
cat("Second-confirmation common-random-number grid study completed.\n")
