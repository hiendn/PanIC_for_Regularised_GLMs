## Confirmatory PanIC-CF calibration and comparison functions.
##
## This file sources the locked preceding implementation and replaces the
## undershoot weight and replication driver.  The assessed methods are
## PanIC-CF, ordinary minimum-loss five-fold CV, and a family-appropriate
## BIC-like comparator.

candidate_source_dir <- local({
  candidates <- vapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.null(value)) NA_character_ else as.character(value)
  }, character(1))
  candidates <- candidates[!is.na(candidates)]
  if (length(candidates)) {
    dirname(normalizePath(tail(candidates, 1L)))
  } else {
    getwd()
  }
})
prior_production_dir <- normalizePath(
  file.path(candidate_source_dir, "reference_implementation"),
  mustWork = TRUE
)
source(file.path(prior_production_dir, "revised_cf_functions.R"))

calibration_weight <- function(n_validation, config = CONFIG) {
  n_validation <- as.numeric(n_validation)
  if (length(n_validation) != 1L || !is.finite(n_validation) ||
      n_validation <= 0) {
    stop("n_validation must be one positive finite value")
  }
  base_weight <- log(log(n_validation + exp(exp(1))))
  base_weight^config$calibration_weight_power
}

## glmnet returns exact zeros for inactive coefficients.  Interpolated path
## points are classified by their literal fitted coefficient values, without
## a post hoc numerical threshold.
support_metrics <- function(beta_hat, truth) {
  selected <- beta_hat != 0
  fp <- sum(selected & !truth$support)
  fn <- sum(!selected & truth$support)
  list(
    fp = fp,
    fn = fn,
    fpr = fp / sum(!truth$support),
    fnr = fn / sum(truth$support),
    exact = as.integer(fp == 0L && fn == 0L),
    wrong = fp + fn,
    norm_error = sum(abs(beta_hat)) - truth$true_radius
  )
}

evaluate_revised_index <- function(path, index, method, truth, test, family,
                                   kappa = NA_real_,
                                   calibration_failed = FALSE,
                                   default_used = FALSE,
                                   config = CONFIG) {
  beta_hat <- path$beta[, index]
  beta0_hat <- path$beta0[index]
  support <- support_metrics(beta_hat, truth)
  eta_test <- drop(beta0_hat + test$x %*% beta_hat)
  data.frame(
    method = method,
    fp = support$fp,
    fn = support$fn,
    fpr = support$fpr,
    fnr = support$fnr,
    exact = support$exact,
    wrong = support$wrong,
    signed_attained_radius_error = support$norm_error,
    test_deviance = mean_deviance(family, test$y, eta_test),
    selected_grid_radius = path$grid_radius[index],
    attained_radius = path$attained_radius[index],
    selected_lower_endpoint = as.integer(index == 1L),
    selected_upper_endpoint = as.integer(index == length(path$grid_radius)),
    kappa = kappa,
    method_failed = 0L,
    failure_reason = "",
    calibration_failed = as.integer(calibration_failed),
    default_used = as.integer(default_used),
    stringsAsFactors = FALSE
  )
}

calibration_failure_rows <- function(common, split_seeds, n) {
  rows <- vector("list", 2L * length(split_seeds))
  counter <- 1L
  for (repeat_index in seq_along(split_seeds)) {
    for (direction in 1:2) {
      rows[[counter]] <- cbind(
        common,
        data.frame(
          split_repeat = repeat_index, direction = direction,
          split_seed = split_seeds[repeat_index], n_training = n %/% 2L,
          n_validation = n - n %/% 2L, training_path_ok = 0L,
          training_pilot_ok = 0L, validation_pilot_ok = 0L,
          raw_cross_signed_radius = NA_real_,
          projected_cross_signed_radius = NA_real_,
          projected_below_grid = NA_integer_,
          projected_above_grid = NA_integer_, projection_applied = NA_integer_,
          weight = NA_real_,
          training_path_warning_count = NA_integer_,
          training_pilot_warning_count = NA_integer_,
          validation_pilot_warning_count = NA_integer_,
          max_active_radius_error = NA_real_,
          failure_reason =
            "not evaluated because the full-sample path failed",
          stringsAsFactors = FALSE
        )
      )
      counter <- counter + 1L
    }
  }
  do.call(rbind, rows)
}

select_cv <- function(validation_loss, radii) {
  if (!is.matrix(validation_loss) || ncol(validation_loss) != length(radii)) {
    stop("CV loss matrix and radius grid are incompatible")
  }
  if (nrow(validation_loss) < 2L || any(!is.finite(validation_loss))) {
    stop("CV selection requires at least two complete finite fold-loss rows")
  }
  if (is.unsorted(radii, strictly = TRUE)) {
    stop("The radius grid must be strictly increasing")
  }
  mean_loss_by_radius <- colMeans(validation_loss)
  selected_index <- smallest_minimiser(mean_loss_by_radius)
  list(
    index = selected_index,
    mean_at_selected = mean_loss_by_radius[selected_index],
    mean_loss = mean_loss_by_radius
  )
}

unexpected_candidate_error_result <- function(scenario, scenario_index,
                                              replication_index, message,
                                              config = CONFIG) {
  streams <- seed_streams(scenario_index, replication_index)
  truth <- scenario_truth(
    scenario$family[[1L]], scenario$rho[[1L]], config
  )
  common <- make_common_row(scenario, replication_index, streams, truth)
  family <- scenario$family[[1L]]
  bic_label <- "BIC-like"
  reason <- paste0("unexpected error: ", message)
  primary <- rbind(
    empty_evaluation_row("PanIC-CF", TRUE, reason, TRUE, FALSE),
    empty_evaluation_row(bic_label, TRUE, reason, FALSE, FALSE),
    empty_evaluation_row(config$primary_cv_method, TRUE, reason, FALSE, FALSE)
  )
  primary <- cbind(common[rep(1L, nrow(primary)), ], primary)
  calibration <- calibration_failure_rows(
    common, streams$calibration_splits, as.integer(scenario$n[[1L]])
  )
  calibration$failure_reason <- reason
  diagnostic <- cbind(
    common,
    data.frame(
      full_path_failed = 1L, calibration_failed = 1L, cv_failed = 1L,
      full_path_failure_reason = reason,
      calibration_failure_reason = reason, cv_failure_reason = reason,
      calibration_rows_expected = 2L * config$calibration_half_splits,
      calibration_rows_successful = 0L,
      raw_target_below_grid = NA_integer_, raw_target_above_grid = NA_integer_,
      projection_count = NA_integer_, raw_target_min = NA_real_,
      raw_target_mean = NA_real_, raw_target_max = NA_real_,
      projected_target_mean = NA_real_, kappa_hat = NA_real_,
      kappa_lower_endpoint = NA_integer_, kappa_upper_endpoint = NA_integer_,
      default_used = 0L,
      calibration_mean_min = NA_real_, calibration_se_at_min = NA_real_,
      eligible_kappa_count = NA_integer_, unregularised_radius = NA_real_,
      cv_index = NA_integer_, cv_mean_loss_at_selected = NA_real_,
      full_warning_count = NA_integer_, calibration_warning_count = NA_integer_,
      cv_warning_count = NA_integer_, maximum_active_radius_error = NA_real_,
      elapsed_seconds = NA_real_, stringsAsFactors = FALSE
    )
  )
  list(primary = primary, diagnostic = diagnostic, calibration = calibration)
}

run_candidate_replication <- function(scenario, scenario_index,
                                      replication_index,
                                      radii = radius_grid(), config = CONFIG) {
  family <- scenario$family[[1L]]
  n <- as.integer(scenario$n[[1L]])
  rho <- scenario$rho[[1L]]
  streams <- seed_streams(scenario_index, replication_index)
  truth <- scenario_truth(family, rho, config)
  common <- make_common_row(scenario, replication_index, streams, truth)
  train <- generate_dataset(n, family, truth, streams$training, config)
  test <- generate_dataset(config$n_test, family, truth, streams$test, config)
  start_time <- proc.time()[["elapsed"]]

  full_path <- fit_radius_path(train$x, train$y, family, radii, config)
  if (!full_path$ok) {
    reason <- paste0("full-sample path: ", full_path$reason)
    method_names <- c("PanIC-CF", "BIC-like", config$primary_cv_method)
    primary <- do.call(rbind, lapply(method_names, function(method) {
      empty_evaluation_row(
        method, TRUE, reason, method == "PanIC-CF", FALSE
      )
    }))
    calibration_rows <- calibration_failure_rows(
      common, streams$calibration_splits, n
    )
    diagnostics <- cbind(
      common,
      data.frame(
        full_path_failed = 1L, calibration_failed = 1L, cv_failed = 1L,
        full_path_failure_reason = reason,
        calibration_failure_reason = reason, cv_failure_reason = reason,
        calibration_rows_expected =
          2L * config$calibration_half_splits,
        calibration_rows_successful = 0L,
        raw_target_below_grid = NA_integer_,
        raw_target_above_grid = NA_integer_, projection_count = NA_integer_,
        raw_target_min = NA_real_, raw_target_mean = NA_real_,
        raw_target_max = NA_real_, projected_target_mean = NA_real_,
        kappa_hat = NA_real_, kappa_lower_endpoint = NA_integer_,
        kappa_upper_endpoint = NA_integer_,
        default_used = 0L,
        calibration_mean_min = NA_real_, calibration_se_at_min = NA_real_,
        eligible_kappa_count = NA_integer_, unregularised_radius = NA_real_,
        cv_index = NA_integer_, cv_mean_loss_at_selected = NA_real_,
        full_warning_count = length(full_path$warnings),
        calibration_warning_count = NA_integer_,
        cv_warning_count = NA_integer_,
        maximum_active_radius_error = NA_real_,
        elapsed_seconds = proc.time()[["elapsed"]] - start_time,
        stringsAsFactors = FALSE
      )
    )
    primary <- cbind(common[rep(1L, nrow(primary)), ], primary)
    return(list(primary = primary, diagnostic = diagnostics,
                calibration = calibration_rows))
  }

  ## The five balanced half-splits, two scoring directions, cross-signed
  ## target, projection, multiplier grid, and one-SE multiplier rule are
  ## unchanged from the preceding locked implementation.
  n_units <- config$calibration_half_splits * config$calibration_directions
  selected_psi <- matrix(
    NA_real_, nrow = n_units, ncol = length(config$kappa_grid)
  )
  discrepancy <- matrix(
    NA_real_, nrow = n_units, ncol = length(config$kappa_grid)
  )
  calibration_rows <- vector("list", n_units)
  calibration_reasons <- character()
  calibration_warning_count <- 0L
  psi <- (radii - min(radii)) / (max(radii) - min(radii))
  unit <- 1L

  for (repeat_index in seq_len(config$calibration_half_splits)) {
    assignment <- balanced_folds(
      n, 2L, streams$calibration_splits[repeat_index]
    )
    for (direction in 1:2) {
      validation_label <- direction
      validation <- which(assignment == validation_label)
      training <- which(assignment != validation_label)
      fold_path <- fit_radius_path(
        train$x[training, , drop = FALSE], train$y[training], family,
        radii, config
      )
      training_pilot <- fit_unregularised(
        train$x[training, , drop = FALSE], train$y[training], family, config
      )
      validation_pilot <- fit_unregularised(
        train$x[validation, , drop = FALSE], train$y[validation], family,
        config
      )
      reasons <- character()
      if (!fold_path$ok) {
        reasons <- c(reasons, paste0("training path: ", fold_path$reason))
      }
      if (!training_pilot$ok) {
        reasons <- c(
          reasons, paste0("training pilot: ", training_pilot$reason)
        )
      }
      if (!validation_pilot$ok) {
        reasons <- c(
          reasons, paste0("validation pilot: ", validation_pilot$reason)
        )
      }
      warning_count <- length(fold_path$warnings) +
        length(training_pilot$warnings) + length(validation_pilot$warnings)
      calibration_warning_count <- calibration_warning_count + warning_count
      raw_target <- projected_target <- NA_real_
      below <- above <- projection <- NA_integer_
      weight <- calibration_weight(length(validation), config)

      if (!length(reasons)) {
        raw_target <- cross_signed_radius(
          training_pilot$beta, validation_pilot$beta
        )
        projected_target <- project_radius(
          raw_target, min(radii), max(radii)
        )
        below <- as.integer(raw_target < min(radii))
        above <- as.integer(raw_target > max(radii))
        projection <- as.integer(below == 1L || above == 1L)
        target_psi <- (projected_target - min(radii)) /
          (max(radii) - min(radii))
        penalty_shape <- (1 + psi) *
          sqrt(log(length(training)) / length(training))
        for (ell in seq_along(config$kappa_grid)) {
          path_index <- smallest_minimiser(
            fold_path$risk + config$kappa_grid[ell] * penalty_shape
          )
          selected_psi[unit, ell] <- psi[path_index]
          discrepancy[unit, ell] <-
            pmax(selected_psi[unit, ell] - target_psi, 0)^2 +
            weight * pmax(target_psi - selected_psi[unit, ell], 0)^2
        }
      } else {
        calibration_reasons <- c(
          calibration_reasons,
          paste0(
            "split ", repeat_index, ", direction ", direction, ": ",
            paste(reasons, collapse = "; ")
          )
        )
      }

      calibration_rows[[unit]] <- cbind(
        common,
        data.frame(
          split_repeat = repeat_index, direction = direction,
          split_seed = streams$calibration_splits[repeat_index],
          n_training = length(training), n_validation = length(validation),
          training_path_ok = as.integer(fold_path$ok),
          training_pilot_ok = as.integer(training_pilot$ok),
          validation_pilot_ok = as.integer(validation_pilot$ok),
          raw_cross_signed_radius = raw_target,
          projected_cross_signed_radius = projected_target,
          projected_below_grid = below, projected_above_grid = above,
          projection_applied = projection, weight = weight,
          training_path_warning_count = length(fold_path$warnings),
          training_pilot_warning_count = length(training_pilot$warnings),
          validation_pilot_warning_count = length(validation_pilot$warnings),
          max_active_radius_error = if (fold_path$ok) {
            fold_path$max_active_radius_error
          } else {
            NA_real_
          },
          failure_reason = paste(reasons, collapse = "; "),
          stringsAsFactors = FALSE
        )
      )
      unit <- unit + 1L
    }
  }

  calibration_rows <- do.call(rbind, calibration_rows)
  calibration_failed <- length(calibration_reasons) > 0L
  default_used <- FALSE
  if (!calibration_failed) {
    selection <- select_largest_kappa_one_se(
      discrepancy, config$kappa_grid
    )
  } else {
    if (min(abs(
      config$kappa_grid - config$calibration_default_kappa
    )) > 1e-12) {
      stop("The prespecified calibration default is not on the kappa grid")
    }
    default_index <- which.min(abs(
      config$kappa_grid - config$calibration_default_kappa
    ))
    selection <- list(
      index = default_index,
      kappa = config$kappa_grid[default_index], min_index = NA_integer_,
      mean_at_min = NA_real_, se_at_min = NA_real_,
      eligible_count = NA_integer_,
      lower_endpoint = as.integer(default_index == 1L),
      upper_endpoint = as.integer(
        default_index == length(config$kappa_grid)
      )
    )
    default_used <- TRUE
  }
  full_penalty_shape <- (1 + psi) * sqrt(log(n) / n)
  pan_index <- smallest_minimiser(
    full_path$risk + selection$kappa * full_penalty_shape
  )
  ## Ordinary five-fold CV selects the smallest radius attaining the minimum
  ## mean validation loss on the common radius grid.
  cv_assignment <- balanced_folds(n, config$cv_folds, streams$cv_folds)
  cv_validation_loss <- matrix(
    NA_real_, nrow = config$cv_folds, ncol = length(radii)
  )
  cv_reasons <- character()
  cv_warning_count <- 0L
  cv_radius_error <- rep(NA_real_, config$cv_folds)
  for (fold in seq_len(config$cv_folds)) {
    validation <- which(cv_assignment == fold)
    training <- which(cv_assignment != fold)
    fold_path <- fit_radius_path(
      train$x[training, , drop = FALSE], train$y[training], family,
      radii, config
    )
    cv_warning_count <- cv_warning_count + length(fold_path$warnings)
    if (!fold_path$ok) {
      cv_reasons <- c(
        cv_reasons, paste0("fold ", fold, ": ", fold_path$reason)
      )
    } else {
      eta_validation <- sweep(
        train$x[validation, , drop = FALSE] %*% fold_path$beta,
        2L, fold_path$beta0, FUN = "+"
      )
      cv_validation_loss[fold, ] <- mean_loss(
        family, train$y[validation], eta_validation
      )
      cv_radius_error[fold] <- fold_path$max_active_radius_error
    }
  }
  cv_failed <- length(cv_reasons) > 0L
  cv_selection <- if (cv_failed) {
    NULL
  } else {
    select_cv(cv_validation_loss, radii)
  }

  ## Support is the literal nonzero pattern returned by the fitted path.  The
  ## BIC-like comparator uses the same exact active-set definition.
  active_count <- colSums(full_path$beta != 0)
  monotone_active_count <- cummax(active_count)
  if (family == "gaussian") {
    bic_criterion <- full_path$risk +
      (monotone_active_count + config$bic_epsilon * psi) * log(n) / n
    bic_label <- "BIC-like"
  } else {
    bic_criterion <- full_path$risk +
      0.5 * (monotone_active_count + config$bic_epsilon * psi) *
      log(n) / n
    bic_label <- "BIC-like"
  }
  bic_index <- smallest_minimiser(bic_criterion)

  cv_failure_row <- function(method) {
    empty_evaluation_row(
      method, TRUE, paste(cv_reasons, collapse = " | "), FALSE, FALSE
    )
  }
  primary <- rbind(
    evaluate_revised_index(
      full_path, pan_index, "PanIC-CF", truth, test, family,
      selection$kappa, calibration_failed, default_used, config
    ),
    evaluate_revised_index(
      full_path, bic_index, bic_label, truth, test, family,
      NA_real_, FALSE, FALSE, config
    ),
    if (cv_failed) {
      cv_failure_row(config$primary_cv_method)
    } else {
      evaluate_revised_index(
        full_path, cv_selection$index, config$primary_cv_method,
        truth, test, family, NA_real_, FALSE, FALSE, config
      )
    }
  )
  primary <- cbind(common[rep(1L, nrow(primary)), ], primary)

  raw_target <- calibration_rows$raw_cross_signed_radius
  projected_target <- calibration_rows$projected_cross_signed_radius
  successful_rows <- is.finite(raw_target)
  maximum_radius_error <- max(
    c(
      full_path$max_active_radius_error,
      calibration_rows$max_active_radius_error,
      cv_radius_error
    ),
    na.rm = TRUE
  )
  diagnostics <- cbind(
    common,
    data.frame(
      full_path_failed = 0L,
      calibration_failed = as.integer(calibration_failed),
      cv_failed = as.integer(cv_failed), full_path_failure_reason = "",
      calibration_failure_reason = paste(
        calibration_reasons, collapse = " | "
      ),
      cv_failure_reason = paste(cv_reasons, collapse = " | "),
      calibration_rows_expected = n_units,
      calibration_rows_successful = sum(successful_rows),
      raw_target_below_grid = sum(raw_target < min(radii), na.rm = TRUE),
      raw_target_above_grid = sum(raw_target > max(radii), na.rm = TRUE),
      projection_count = sum(
        calibration_rows$projection_applied, na.rm = TRUE
      ),
      raw_target_min = if (any(successful_rows)) {
        min(raw_target[successful_rows])
      } else {
        NA_real_
      },
      raw_target_mean = if (any(successful_rows)) {
        mean(raw_target[successful_rows])
      } else {
        NA_real_
      },
      raw_target_max = if (any(successful_rows)) {
        max(raw_target[successful_rows])
      } else {
        NA_real_
      },
      projected_target_mean = if (any(successful_rows)) {
        mean(projected_target[successful_rows])
      } else {
        NA_real_
      },
      kappa_hat = selection$kappa,
      kappa_lower_endpoint = selection$lower_endpoint,
      kappa_upper_endpoint = selection$upper_endpoint,
      default_used = as.integer(default_used),
      calibration_mean_min = selection$mean_at_min,
      calibration_se_at_min = selection$se_at_min,
      eligible_kappa_count = selection$eligible_count,
      unregularised_radius = full_path$unregularised_radius,
      cv_index = if (cv_failed) NA_integer_ else cv_selection$index,
      cv_mean_loss_at_selected = if (cv_failed) {
        NA_real_
      } else {
        cv_selection$mean_at_selected
      },
      full_warning_count = length(full_path$warnings),
      calibration_warning_count = calibration_warning_count,
      cv_warning_count = cv_warning_count,
      maximum_active_radius_error = maximum_radius_error,
      elapsed_seconds = proc.time()[["elapsed"]] - start_time,
      stringsAsFactors = FALSE
    )
  )
  list(primary = primary, diagnostic = diagnostics,
       calibration = calibration_rows)
}

rm(candidate_source_dir, prior_production_dir)
