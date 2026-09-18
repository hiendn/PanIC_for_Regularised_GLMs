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
} else normalizePath("Grid_Sensitivity_Study.R")
production_dir <- dirname(script_path)
source(file.path(production_dir, "Simulation_Config.R"))

## Prespecified, disjoint stream for the common-random-number grid study.
CONFIRMATORY_MASTER_SEED <- CONFIG$master_seed
SENSITIVITY_MASTER_SEED <- 2036091802L
CONFIG$master_seed <- SENSITIVITY_MASTER_SEED
source(file.path(production_dir, "PanIC_CF_Functions.R"))

n_rep <- as.integer(arg_value("--n-rep", 500L))
cores <- max(1L, as.integer(arg_value("--cores", CONFIG$requested_cores)))
output_name <- arg_value("--output-dir", "results")
results_dir <- if (grepl("^/", output_name)) output_name else
  file.path(production_dir, output_name)
generated_dir <- results_dir
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

grid_sizes <- c(61L, 121L, 241L)
scenario <- data.frame(
  scenario_id = "linear_iid_n1000_grid_sensitivity",
  family = "gaussian", n = 1000L, rho = 0,
  stringsAsFactors = FALSE
)
scenario_seed_index <- 1L

for (m in grid_sizes) {
  message("Running revised grid sensitivity m=", m,
          " with common random numbers")
  worker <- function(replication) {
    tryCatch(
      run_revised_replication(
        scenario, scenario_seed_index, replication,
        radii = radius_grid(m), config = CONFIG
      ),
      error = function(e) unexpected_error_result(
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
  } else lapply(seq_len(n_rep), worker)
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
    out$kappa_lower_endpoint_selections <- sum(
      d$kappa_lower_endpoint, na.rm = TRUE
    )
    out$kappa_upper_endpoint_selections <- sum(
      d$kappa_upper_endpoint, na.rm = TRUE
    )
    out$mean_shared_runtime <- mean(d$elapsed_seconds, na.rm = TRUE)
    out$mcse_shared_runtime <- mcse(d$elapsed_seconds)
    out$maximum_radius_interpolation_error <- max(
      d$maximum_active_radius_error, na.rm = TRUE
    )
    out
  }
)
grid_summary <- do.call(rbind, summary_rows)
method_order <- c("PanIC-CF", "BIC-like", "CV")
grid_summary <- grid_summary[
  order(grid_summary$radius_points,
        match(grid_summary$method, method_order)),
]
write.csv(grid_summary,
          file.path(results_dir, "grid_sensitivity_summary.csv"),
          row.names = FALSE)

## Paired differences exploit the common random numbers across grid sizes.
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
paired_contrasts <- do.call(rbind, paired_rows)
write.csv(
  paired_contrasts,
  file.path(results_dir, "grid_sensitivity_paired_contrasts.csv"),
  row.names = FALSE
)

## Check that the random data, split streams, and target GLMs are truly common.
calibration_key <- c("replication", "split_repeat", "direction")
reference <- calibration[calibration$radius_points == 121L,
                         c(calibration_key, "training_seed", "split_seed",
                           "test_seed", "raw_cross_signed_radius")]
for (m in c(61L, 241L)) {
  candidate <- calibration[calibration$radius_points == m,
                           c(calibration_key, "training_seed", "split_seed",
                             "test_seed", "raw_cross_signed_radius")]
  comparison <- merge(reference, candidate, by = calibration_key,
                      suffixes = c("_reference", "_candidate"))
  if (nrow(comparison) != 10L * n_rep ||
      any(comparison$training_seed_reference !=
          comparison$training_seed_candidate) ||
      any(comparison$split_seed_reference != comparison$split_seed_candidate) ||
      any(comparison$test_seed_reference != comparison$test_seed_candidate) ||
      max(abs(comparison$raw_cross_signed_radius_reference -
              comparison$raw_cross_signed_radius_candidate), na.rm = TRUE) >
        1e-12) {
    stop("Common-random-number grid comparison failed for m=", m)
  }
}

cell <- function(mean, se, digits = 3L) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean, se)
}
table_lines <- c(
  "\\begin{table}[H]", "\\centering",
  paste0("\\caption{Grid sensitivity for the revised PanIC-CF procedure in ",
         "the independent Gaussian setting with $n=1000$ and common random ",
         "numbers across $N=500$ replications. Every grid spans $[0,20]$. ",
         "Entries are means with Monte Carlo standard errors in parentheses. ",
         "Runtime is the shared end-to-end time for the full-sample path, ten ",
         "half-sample calibration directions, and five ordinary CV folds.}"),
  "\\label{Table: grid sensitivity}",
  "\\begin{adjustbox}{width=\\textwidth}", "\\scriptsize",
  "\\begin{tabular}{rrlrrrrr}", "\\hline",
  paste0("$m$ & Spacing & Method & Total support error & Exact & Test deviance ",
         "& Selected radius & Runtime (s) \\\\"),
  "\\hline"
)
for (i in seq_len(nrow(grid_summary))) {
  row <- grid_summary[i, ]
  table_lines <- c(table_lines, paste0(
    row$radius_points, " & ", sprintf("%.3f", row$grid_spacing), " & ",
    row$method, " & ", cell(row$wrong_mean, row$wrong_mcse, 2), " & ",
    cell(row$exact_mean, row$exact_mcse, 3), " & ",
    cell(row$test_deviance_mean, row$test_deviance_mcse, 3), " & ",
    cell(row$selected_grid_radius_mean,
         row$selected_grid_radius_mcse, 3), " & ",
    cell(row$mean_shared_runtime, row$mcse_shared_runtime, 3), " \\\\"
  ))
}
table_lines <- c(table_lines, "\\hline", "\\end{tabular}",
                 "\\end{adjustbox}", "\\end{table}")
writeLines(
  table_lines,
  file.path(generated_dir, "table_grid_sensitivity.tex")
)

write.csv(
  data.frame(
    item = c(
      "sensitivity_master_seed", "confirmatory_master_seed",
      "scenario_seed_index", "replications", "workers",
      "radius_points", "radius_interval"
    ),
    value = c(
      SENSITIVITY_MASTER_SEED, CONFIRMATORY_MASTER_SEED,
      scenario_seed_index, n_rep, cores,
      paste(grid_sizes, collapse = ";"),
      paste0("[", CONFIG$radius_min, ",", CONFIG$radius_max, "]")
    ),
    stringsAsFactors = FALSE
  ),
  file.path(results_dir, "grid_sensitivity_configuration.csv"),
  row.names = FALSE
)
cat("Revised common-random-number grid-sensitivity study completed.\n")
