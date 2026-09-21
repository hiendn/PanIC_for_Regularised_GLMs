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
} else normalizePath("benchmark_runtime.R")
production_dir <- dirname(script_path)
.libPaths(c(file.path(production_dir, "Rlib"), .libPaths()))
source(file.path(production_dir, "Simulation_Config.R"))
source(file.path(production_dir, "PanIC_CF_Functions.R"))
## The frozen implementation sourced above defines its own `production_dir`.
## Restore this script's directory before resolving output and renderer paths.
production_dir <- dirname(script_path)
output_name <- arg_value("--output-dir", "results")
results_dir <- if (grepl("^/", output_name)) output_name else
  file.path(production_dir, output_name)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

runtime_seed <- function(family_index, replication, warmup = FALSE) {
  warmup_offset <- if (warmup) 500000L else 0L
  as.integer(CONFIG$runtime_master_seed + warmup_offset +
               100000L * family_index + abs(as.integer(replication)))
}

runtime_dataset <- function(family, family_index, replication,
                            warmup = FALSE) {
  seed <- runtime_seed(family_index, replication, warmup)
  truth <- scenario_truth(family, 0, CONFIG)
  data <- generate_dataset(CONFIG$runtime_n, family, truth, seed, CONFIG)
  list(data = data, truth = truth, seed = seed)
}

time_bic <- function(data, family, radii) {
  psi <- (radii - min(radii)) / (max(radii) - min(radii))
  result <- NULL
  elapsed <- system.time({
    full_path <- fit_radius_path(data$x, data$y, family, radii, CONFIG)
    if (!full_path$ok) stop("full path: ", full_path$reason)
    active_count <- cummax(
      colSums(full_path$beta != 0)
    )
    scale <- if (family == "gaussian") 1 else 0.5
    criterion <- full_path$risk + scale *
      (active_count + CONFIG$bic_epsilon * psi) *
      log(CONFIG$runtime_n) / CONFIG$runtime_n
    result <- smallest_minimiser(criterion)
  })[["elapsed"]]
  list(elapsed = elapsed, selected_index = result,
       paths_fitted = 1L, unregularised_sign_fits = 0L,
       unregularised_validation_fits = 0L)
}

time_cv <- function(data, family, radii, fold_seed) {
  result <- NULL
  elapsed <- system.time({
    full_path <- fit_radius_path(data$x, data$y, family, radii, CONFIG)
    if (!full_path$ok) stop("full path: ", full_path$reason)
    assignment <- balanced_folds(
      CONFIG$runtime_n, CONFIG$cv_folds, fold_seed
    )
    validation_loss <- matrix(
      NA_real_, CONFIG$cv_folds, length(radii)
    )
    for (fold in seq_len(CONFIG$cv_folds)) {
      validation <- which(assignment == fold)
      training <- which(assignment != fold)
      path <- fit_radius_path(
        data$x[training, , drop = FALSE], data$y[training], family,
        radii, CONFIG
      )
      if (!path$ok) stop("CV fold ", fold, ": ", path$reason)
      eta <- sweep(
        data$x[validation, , drop = FALSE] %*% path$beta,
        2L, path$beta0, FUN = "+"
      )
      validation_loss[fold, ] <- mean_loss(
        family, data$y[validation], eta
      )
    }
    result <- smallest_minimiser(colMeans(validation_loss))
  })[["elapsed"]]
  list(elapsed = elapsed, selected_index = result,
       paths_fitted = 1L + CONFIG$cv_folds,
       unregularised_sign_fits = 0L,
       unregularised_validation_fits = 0L)
}

time_revised_panic_cf <- function(data, family, radii, split_seeds) {
  psi <- (radii - min(radii)) / (max(radii) - min(radii))
  selected_index <- NULL
  selected_kappa <- NULL
  elapsed <- system.time({
    full_path <- fit_radius_path(data$x, data$y, family, radii, CONFIG)
    if (!full_path$ok) stop("full path: ", full_path$reason)
    n_units <- CONFIG$calibration_half_splits *
      CONFIG$calibration_directions
    discrepancy <- matrix(
      NA_real_, n_units, length(CONFIG$kappa_grid)
    )
    unit <- 1L
    for (repeat_index in seq_len(CONFIG$calibration_half_splits)) {
      assignment <- balanced_folds(
        CONFIG$runtime_n, 2L, split_seeds[repeat_index]
      )
      for (direction in 1:2) {
        validation <- which(assignment == direction)
        training <- which(assignment != direction)
        path <- fit_radius_path(
          data$x[training, , drop = FALSE], data$y[training], family,
          radii, CONFIG
        )
        if (!path$ok) {
          stop("half split ", repeat_index, " direction ", direction,
               " path: ", path$reason)
        }
        training_pilot <- fit_unregularised(
          data$x[training, , drop = FALSE], data$y[training], family,
          CONFIG
        )
        validation_pilot <- fit_unregularised(
          data$x[validation, , drop = FALSE], data$y[validation], family,
          CONFIG
        )
        if (!training_pilot$ok) {
          stop("training-sign GLM: ", training_pilot$reason)
        }
        if (!validation_pilot$ok) {
          stop("validation GLM: ", validation_pilot$reason)
        }
        raw_target <- cross_signed_radius(
          training_pilot$beta, validation_pilot$beta
        )
        target <- project_radius(raw_target, min(radii), max(radii))
        target_psi <- (target - min(radii)) /
          (max(radii) - min(radii))
        shape <- (1 + psi) * sqrt(log(length(training)) / length(training))
        weight <- calibration_weight(length(validation))
        for (ell in seq_along(CONFIG$kappa_grid)) {
          index <- smallest_minimiser(
            path$risk + CONFIG$kappa_grid[ell] * shape
          )
          selected <- psi[index]
          discrepancy[unit, ell] <-
            pmax(selected - target_psi, 0)^2 +
            weight * pmax(target_psi - selected, 0)^2
        }
        unit <- unit + 1L
      }
    }
    selection <- select_largest_kappa_one_se(
      discrepancy, CONFIG$kappa_grid
    )
    selected_kappa <- selection$kappa
    full_shape <- (1 + psi) *
      sqrt(log(CONFIG$runtime_n) / CONFIG$runtime_n)
    selected_index <- smallest_minimiser(
      full_path$risk + selected_kappa * full_shape
    )
  })[["elapsed"]]
  list(
    elapsed = elapsed, selected_index = selected_index,
    selected_kappa = selected_kappa,
    paths_fitted = 1L +
      CONFIG$calibration_half_splits * CONFIG$calibration_directions,
    unregularised_sign_fits =
      CONFIG$calibration_half_splits * CONFIG$calibration_directions,
    unregularised_validation_fits =
      CONFIG$calibration_half_splits * CONFIG$calibration_directions
  )
}

benchmark_method <- function(family, family_index, method, replication,
                             warmup = FALSE) {
  generated <- runtime_dataset(
    family, family_index, replication, warmup
  )
  data <- generated$data
  seed <- generated$seed
  radii <- radius_grid()
  fit <- tryCatch(
    switch(
      method,
      `PanIC-CF` = time_revised_panic_cf(
        data, family, radii,
        as.integer(seed + 10000L +
          1000L * (seq_len(CONFIG$calibration_half_splits) - 1L))
      ),
      `CV` = time_cv(data, family, radii, seed + 30000L),
      `BIC-like` = time_bic(data, family, radii),
      stop("Unknown method: ", method)
    ),
    error = function(e) list(
      elapsed = NA_real_, selected_index = NA_integer_,
      selected_kappa = NA_real_, paths_fitted = NA_integer_,
      unregularised_sign_fits = NA_integer_,
      unregularised_validation_fits = NA_integer_,
      error = conditionMessage(e)
    )
  )
  data.frame(
    family = family, n = CONFIG$runtime_n, method = method,
    replication = as.integer(replication), data_seed = seed,
    elapsed_seconds = fit$elapsed,
    selected_index = fit$selected_index,
    selected_kappa = if (is.null(fit$selected_kappa)) NA_real_ else
      fit$selected_kappa,
    paths_fitted = fit$paths_fitted,
    unregularised_sign_fits = fit$unregularised_sign_fits,
    unregularised_validation_fits = fit$unregularised_validation_fits,
    failed = as.integer(!is.finite(fit$elapsed)),
    failure_reason = if (is.null(fit$error)) "" else fit$error,
    stringsAsFactors = FALSE
  )
}

families <- c("gaussian", "binomial", "poisson")
methods <- c("PanIC-CF", "CV", "BIC-like")

for (family_index in seq_along(families)) {
  family <- families[family_index]
  message("Warming runtime benchmark for ", family)
  for (warmup in seq_len(CONFIG$runtime_warmup)) {
    for (method in methods) {
      warm <- benchmark_method(
        family, family_index, method, warmup, warmup = TRUE
      )
      if (warm$failed) stop("Runtime warmup failed: ", warm$failure_reason)
    }
  }
}

rows <- vector(
  "list", length(families) * length(methods) *
    CONFIG$runtime_replications
)
cursor <- 1L
for (family_index in seq_along(families)) {
  family <- families[family_index]
  message("Timing ", family)
  for (replication in seq_len(CONFIG$runtime_replications)) {
    for (method in methods) {
      rows[[cursor]] <- benchmark_method(
        family, family_index, method, replication, warmup = FALSE
      )
      cursor <- cursor + 1L
    }
  }
}
runtime <- do.call(rbind, rows)
write.csv(runtime, file.path(results_dir, "runtime_raw.csv"), row.names = FALSE)

summary_rows <- lapply(
  split(runtime, interaction(runtime$family, runtime$method, drop = TRUE)),
  function(dat) {
    elapsed <- dat$elapsed_seconds[is.finite(dat$elapsed_seconds)]
    data.frame(
      family = dat$family[1L], n = dat$n[1L], method = dat$method[1L],
      attempted_replications = nrow(dat), failures = sum(dat$failed),
      mean = if (length(elapsed)) mean(elapsed) else NA_real_,
      se = if (length(elapsed) > 1L) {
        sd(elapsed) / sqrt(length(elapsed))
      } else NA_real_,
      median = if (length(elapsed)) median(elapsed) else NA_real_,
      minimum = if (length(elapsed)) min(elapsed) else NA_real_,
      maximum = if (length(elapsed)) max(elapsed) else NA_real_,
      paths_fitted = unique(dat$paths_fitted[is.finite(dat$paths_fitted)])[1L],
      unregularised_sign_fits = unique(
        dat$unregularised_sign_fits[
          is.finite(dat$unregularised_sign_fits)
        ]
      )[1L],
      unregularised_validation_fits = unique(
        dat$unregularised_validation_fits[
          is.finite(dat$unregularised_validation_fits)
        ]
      )[1L],
      stringsAsFactors = FALSE
    )
  }
)
runtime_summary <- do.call(rbind, summary_rows)
runtime_summary <- runtime_summary[
  order(match(runtime_summary$family, families),
        match(runtime_summary$method, methods)),
]
write.csv(runtime_summary, file.path(results_dir, "runtime_summary.csv"),
          row.names = FALSE)

source(file.path(production_dir, "Render_Runtime_Figure.R"), local = TRUE)

hardware_raw <- tryCatch(
  system2("system_profiler", "SPHardwareDataType", stdout = TRUE,
          stderr = FALSE),
  error = function(e) character()
)
hardware_summary <- trimws(hardware_raw[grepl(
  "Model Name:|Model Identifier:|Chip:|Total Number of Cores:|Memory:",
  hardware_raw
)])
if (!length(hardware_summary)) hardware_summary <- "Hardware detail unavailable"
writeLines(
  c(
    paste0("Timing date: ", format(Sys.time(),
      tz = "Australia/Melbourne", usetz = TRUE)),
    paste0("Operating system: ", Sys.info()[["sysname"]], " ",
           Sys.info()[["release"]]),
    paste0("Machine: ", Sys.info()[["machine"]]),
    hardware_summary,
    "Timing worker count: 1",
    paste0("Runtime master seed: ", CONFIG$runtime_master_seed),
    paste0("Warm-up runs per family and method: ", CONFIG$runtime_warmup),
    paste0("Timed replications per family and method: ",
           CONFIG$runtime_replications),
    "Elapsed time excludes data generation and independent test evaluation.",
    paste0("PanIC-CF includes ten half-sample path fits, ten ",
           "training-sign GLMs, ten independent validation GLMs, and one ",
           "full-sample path."),
    "CV includes five 80%-sample paths and one full-sample path.",
    paste0("The BIC-like comparator includes one full-sample path and ",
           "family-appropriate active-count criterion evaluation."),
    paste0("R: ", R.version.string),
    paste0("glmnet: ", as.character(packageVersion("glmnet")))
  ),
  file.path(results_dir, "runtime_environment.txt")
)

if (any(runtime$failed)) {
  stop("At least one timed runtime replication failed; see runtime_raw.csv")
}
message("Runtime benchmark and figure completed")
