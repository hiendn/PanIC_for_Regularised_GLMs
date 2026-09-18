## Core simulation and fitting functions.
##
## Candidate models are fixed l1-radius balls.  glmnet is used to compute a
## dense, warm-started penalised path; neighbouring path fits are interpolated
## onto the prespecified radius grid.  A separate validation script compares
## this interpolation with direct radius-targeted root searches.

suppressPackageStartupMessages(library(glmnet))

softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

ar1_covariance <- function(d, rho) {
  outer(seq_len(d), seq_len(d), function(i, j) rho^abs(i - j))
}

scenario_truth <- function(family, rho, config = CONFIG) {
  sigma_x <- ar1_covariance(config$d, rho)
  beta_raw <- numeric(config$d)
  beta_raw[config$active] <- rep(c(1, -1), length.out = length(config$active))
  raw_variance <- drop(crossprod(beta_raw, sigma_x %*% beta_raw))
  beta <- beta_raw * sqrt(config$target_linear_predictor_variance / raw_variance)
  signal_variance <- drop(crossprod(beta, sigma_x %*% beta))
  beta0 <- switch(
    family,
    gaussian = 0,
    binomial = 0,
    poisson = log(config$poisson_marginal_mean) - signal_variance / 2,
    stop("Unsupported family: ", family)
  )
  list(
    beta0 = beta0,
    beta = beta,
    support = abs(beta) > 0,
    sigma_x = sigma_x,
    signal_variance = signal_variance,
    true_radius = sum(abs(beta))
  )
}

generate_dataset <- function(n, family, truth, seed, config = CONFIG) {
  set.seed(seed)
  z <- matrix(rnorm(n * config$d), nrow = n, ncol = config$d)
  x <- z %*% chol(truth$sigma_x)
  eta <- drop(truth$beta0 + x %*% truth$beta)
  y <- switch(
    family,
    gaussian = eta + config$gaussian_sigma * rnorm(n),
    binomial = rbinom(n, size = 1L, prob = plogis(eta)),
    poisson = rpois(n, lambda = exp(eta)),
    stop("Unsupported family: ", family)
  )
  list(x = x, y = y)
}

loss_matrix <- function(family, y, eta) {
  if (is.null(dim(eta))) eta <- matrix(eta, ncol = 1L)
  yy <- matrix(y, nrow = length(y), ncol = ncol(eta))
  switch(
    family,
    gaussian = (yy - eta)^2 / CONFIG$gaussian_sigma^2,
    binomial = softplus(eta) - yy * eta,
    poisson = exp(eta) - yy * eta,
    stop("Unsupported family: ", family)
  )
}

mean_loss <- function(family, y, eta) {
  colMeans(loss_matrix(family, y, eta))
}

mean_deviance <- function(family, y, eta) {
  if (family == "gaussian") {
    return(mean((y - eta)^2 / CONFIG$gaussian_sigma^2))
  }
  if (family == "binomial") {
    return(mean(2 * (softplus(eta) - y * eta)))
  }
  mu <- exp(eta)
  log_term <- numeric(length(y))
  positive <- y > 0
  log_term[positive] <- y[positive] * log(y[positive] / mu[positive])
  mean(2 * (log_term - (y - mu)))
}

intercept_only_fit <- function(y, family) {
  if (family == "gaussian") return(mean(y))
  if (family == "binomial") {
    p <- mean(y)
    if (!is.finite(p) || p <= 0 || p >= 1) return(NA_real_)
    return(qlogis(p))
  }
  mu <- mean(y)
  if (!is.finite(mu) || mu <= 0) return(NA_real_)
  log(mu)
}

fit_unregularised <- function(x, y, family, config = CONFIG) {
  p <- ncol(x) + 1L
  warnings_seen <- character()
  result <- tryCatch(
    withCallingHandlers({
      design <- cbind(`(Intercept)` = 1, x)
      if (family == "gaussian") {
        fit <- lm.fit(design, y, singular.ok = TRUE)
        list(
          coefficients = fit$coefficients,
          rank = fit$rank,
          converged = TRUE,
          boundary = FALSE
        )
      } else {
        fam <- if (family == "binomial") binomial() else poisson()
        fit <- glm.fit(
          x = design, y = y, family = fam, intercept = FALSE,
          control = glm.control(epsilon = 1e-10, maxit = 100L, trace = FALSE)
        )
        list(
          coefficients = fit$coefficients,
          rank = fit$rank,
          converged = isTRUE(fit$converged),
          boundary = isTRUE(fit$boundary)
        )
      }
    }, warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(message = conditionMessage(e)), class = "pilot_error")
  )

  if (inherits(result, "pilot_error")) {
    return(list(ok = FALSE, reason = paste0("error: ", result$message),
                warnings = warnings_seen))
  }
  coef <- as.numeric(result$coefficients)
  warning_failure <- any(grepl(
    "did not converge|fitted probabilities numerically 0 or 1|algorithm stopped",
    warnings_seen, ignore.case = TRUE
  ))
  reasons <- character()
  if (result$rank < p) reasons <- c(reasons, "rank deficient")
  if (!result$converged) reasons <- c(reasons, "did not converge")
  if (result$boundary) reasons <- c(reasons, "boundary fit")
  if (warning_failure) reasons <- c(reasons, "diagnostic warning")
  if (length(coef) != p || any(!is.finite(coef))) reasons <- c(reasons, "non-finite coefficient")
  if (length(coef) == p && any(abs(coef) > config$pilot_coefficient_cap)) {
    reasons <- c(reasons, "coefficient cap exceeded")
  }
  if (length(reasons)) {
    return(list(ok = FALSE, reason = paste(unique(reasons), collapse = "; "),
                warnings = warnings_seen))
  }
  list(
    ok = TRUE,
    beta0 = coef[1L],
    beta = coef[-1L],
    radius = sum(abs(coef[-1L])),
    reason = "",
    warnings = warnings_seen
  )
}

glmnet_path_raw <- function(x, y, family, config = CONFIG) {
  warnings_seen <- character()
  fit <- tryCatch(
    withCallingHandlers(
      glmnet(
        x = x, y = y, family = family, alpha = 1,
        intercept = TRUE, standardize = FALSE,
        nlambda = config$glmnet_nlambda,
        lambda.min.ratio = config$glmnet_lambda_min_ratio,
        thresh = config$glmnet_threshold,
        maxit = config$glmnet_maxit
      ),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(list(message = conditionMessage(e)), class = "path_error")
  )
  if (inherits(fit, "path_error")) {
    return(list(ok = FALSE, reason = paste0("glmnet error: ", fit$message),
                warnings = warnings_seen))
  }
  if (fit$jerr != 0) {
    return(list(ok = FALSE, reason = paste0("glmnet jerr=", fit$jerr),
                warnings = warnings_seen))
  }
  list(ok = TRUE, fit = fit, warnings = warnings_seen)
}

segment_at_radius <- function(beta_lo, beta_hi, target, tol) {
  f <- function(t) sum(abs(beta_lo + t * (beta_hi - beta_lo))) - target
  flo <- f(0)
  fhi <- f(1)
  if (abs(flo) <= tol) return(0)
  if (abs(fhi) <= tol) return(1)
  if (flo > 0 || fhi < 0) {
    stop("Radius target is not bracketed by interpolation segment")
  }
  uniroot(f, interval = c(0, 1), tol = tol)$root
}

fit_radius_path <- function(x, y, family, radii = radius_grid(), config = CONFIG) {
  raw <- glmnet_path_raw(x, y, family, config)
  if (!raw$ok) return(raw)
  fit <- raw$fit
  beta_path <- as.matrix(fit$beta)
  beta0_path <- as.numeric(fit$a0)
  path_radius <- colSums(abs(beta_path))

  b0_zero <- intercept_only_fit(y, family)
  if (!is.finite(b0_zero)) {
    return(list(ok = FALSE, reason = "intercept-only fit does not exist",
                warnings = raw$warnings))
  }
  unreg <- fit_unregularised(x, y, family, config)
  if (!unreg$ok) {
    return(list(ok = FALSE, reason = paste0("full unregularised fit: ", unreg$reason),
                warnings = c(raw$warnings, unreg$warnings)))
  }

  beta_path <- cbind(numeric(ncol(x)), beta_path, unreg$beta)
  beta0_path <- c(b0_zero, beta0_path, unreg$beta0)
  path_radius <- c(0, path_radius, unreg$radius)

  ## Radius is monotone along an exact lasso path.  Sort and remove numerical
  ## duplicates so that interpolation is stable even under machine-level drift.
  ord <- order(path_radius, decreasing = FALSE)
  beta_path <- beta_path[, ord, drop = FALSE]
  beta0_path <- beta0_path[ord]
  path_radius <- path_radius[ord]
  keep <- !duplicated(signif(path_radius, 13L))
  beta_path <- beta_path[, keep, drop = FALSE]
  beta0_path <- beta0_path[keep]
  path_radius <- path_radius[keep]

  m <- length(radii)
  beta <- matrix(0, nrow = ncol(x), ncol = m)
  beta0 <- numeric(m)
  attained <- numeric(m)
  inactive <- logical(m)

  for (j in seq_along(radii)) {
    target <- radii[j]
    if (target <= config$interpolation_radius_tolerance) {
      beta0[j] <- b0_zero
      next
    }
    if (target >= unreg$radius - config$interpolation_radius_tolerance) {
      beta[, j] <- unreg$beta
      beta0[j] <- unreg$beta0
      attained[j] <- unreg$radius
      inactive[j] <- TRUE
      next
    }
    hi <- which(path_radius >= target)[1L]
    if (is.na(hi) || hi <= 1L) {
      return(list(ok = FALSE, reason = paste0("path failed to bracket radius ", target),
                  warnings = raw$warnings))
    }
    lo <- hi - 1L
    tt <- segment_at_radius(
      beta_path[, lo], beta_path[, hi], target,
      config$interpolation_radius_tolerance
    )
    beta[, j] <- beta_path[, lo] + tt * (beta_path[, hi] - beta_path[, lo])
    beta0[j] <- beta0_path[lo] + tt * (beta0_path[hi] - beta0_path[lo])
    attained[j] <- sum(abs(beta[, j]))
  }

  eta <- sweep(x %*% beta, 2L, beta0, FUN = "+")
  risk <- mean_loss(family, y, eta)
  list(
    ok = TRUE,
    beta0 = beta0,
    beta = beta,
    risk = risk,
    grid_radius = radii,
    attained_radius = attained,
    inactive = inactive,
    unregularised_radius = unreg$radius,
    max_active_radius_error = if (any(!inactive)) {
      max(abs(attained[!inactive] - radii[!inactive]))
    } else 0,
    warnings = c(raw$warnings, unreg$warnings),
    glmnet_lambda_count = length(fit$lambda)
  )
}

balanced_folds <- function(n, v, seed) {
  set.seed(seed)
  sample(rep(seq_len(v), length.out = n), size = n, replace = FALSE)
}

smallest_minimiser <- function(x, tolerance = 1e-12) {
  min(which(x <= min(x) + tolerance))
}

support_metrics <- function(beta_hat, truth, active_tolerance = CONFIG$active_tolerance) {
  selected <- abs(beta_hat) > active_tolerance
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

evaluate_path_index <- function(path, index, method, truth, test, family,
                                kappa = NA_real_, calibration_failed = FALSE,
                                config = CONFIG) {
  beta_hat <- path$beta[, index]
  beta0_hat <- path$beta0[index]
  sm <- support_metrics(beta_hat, truth, config$active_tolerance)
  eta_test <- drop(beta0_hat + test$x %*% beta_hat)
  data.frame(
    method = method,
    fp = sm$fp,
    fn = sm$fn,
    fpr = sm$fpr,
    fnr = sm$fnr,
    exact = sm$exact,
    wrong = sm$wrong,
    norm_error = sm$norm_error,
    test_deviance = mean_deviance(family, test$y, eta_test),
    selected_grid_radius = path$grid_radius[index],
    attained_radius = path$attained_radius[index],
    selected_endpoint = as.integer(index == length(path$grid_radius)),
    kappa = kappa,
    calibration_failed = as.integer(calibration_failed),
    stringsAsFactors = FALSE
  )
}

failed_method_row <- function(method) {
  data.frame(
    method = method,
    fp = NA_real_, fn = NA_real_, fpr = NA_real_, fnr = NA_real_,
    exact = NA_real_, wrong = NA_real_, norm_error = NA_real_,
    test_deviance = NA_real_, selected_grid_radius = NA_real_,
    attained_radius = NA_real_, selected_endpoint = NA_real_,
    kappa = NA_real_, calibration_failed = 1L,
    stringsAsFactors = FALSE
  )
}

run_replication <- function(scenario, scenario_index, replication_index,
                            radii = radius_grid(), config = CONFIG) {
  family <- scenario$family[[1L]]
  n <- as.integer(scenario$n[[1L]])
  rho <- scenario$rho[[1L]]
  seed <- replication_seed(scenario_index, replication_index)
  fold_seed <- seed + 30000L
  test_seed <- seed + 60000L
  truth <- scenario_truth(family, rho, config)
  train <- generate_dataset(n, family, truth, seed, config)
  test <- generate_dataset(config$n_test, family, truth, test_seed, config)
  folds <- balanced_folds(n, config$folds, fold_seed)

  fit_start <- proc.time()[["elapsed"]]
  full_path <- fit_radius_path(train$x, train$y, family, radii, config)
  if (!full_path$ok) {
    stop("Full-sample path failed in ", scenario$scenario_id[[1L]],
         " replication ", replication_index, ": ", full_path$reason)
  }

  v <- config$folds
  fold_train_risk <- matrix(NA_real_, nrow = v, ncol = length(radii))
  fold_validation_loss <- matrix(NA_real_, nrow = v, ncol = length(radii))
  pilot_radius <- rep(NA_real_, v)
  pilot_ok <- rep(FALSE, v)
  pilot_reason <- rep("", v)
  fold_warning_count <- integer(v)
  fold_path_radius_error <- numeric(v)

  for (fold in seq_len(v)) {
    validation <- which(folds == fold)
    training <- which(folds != fold)
    fold_path <- fit_radius_path(
      train$x[training, , drop = FALSE], train$y[training], family, radii, config
    )
    if (!fold_path$ok) {
      stop("Fold path failed in ", scenario$scenario_id[[1L]],
           " replication ", replication_index, ", fold ", fold,
           ": ", fold_path$reason)
    }
    fold_train_risk[fold, ] <- fold_path$risk
    eta_validation <- sweep(
      train$x[validation, , drop = FALSE] %*% fold_path$beta,
      2L, fold_path$beta0, FUN = "+"
    )
    fold_validation_loss[fold, ] <- mean_loss(
      family, train$y[validation], eta_validation
    )
    pilot <- fit_unregularised(
      train$x[validation, , drop = FALSE], train$y[validation], family, config
    )
    pilot_ok[fold] <- pilot$ok
    if (pilot$ok) pilot_radius[fold] <- pilot$radius
    if (!pilot$ok) pilot_reason[fold] <- pilot$reason
    fold_warning_count[fold] <- length(c(fold_path$warnings, pilot$warnings))
    fold_path_radius_error[fold] <- fold_path$max_active_radius_error
  }

  psi <- (radii - min(radii)) / (max(radii) - min(radii))
  kappa_grid <- config$kappa_grid
  calibration_failed <- !all(pilot_ok)
  kappa_hat <- NA_real_
  pan_index <- NA_integer_
  kappa_endpoint <- NA_integer_
  calibration_mean_min <- NA_real_
  calibration_se_at_min <- NA_real_

  if (!calibration_failed) {
    discrepancy <- matrix(NA_real_, nrow = v, ncol = length(kappa_grid))
    for (fold in seq_len(v)) {
      n_training <- sum(folds != fold)
      n_validation <- sum(folds == fold)
      penalty_shape <- (1 + psi) * sqrt(log(n_training) / n_training)
      target <- (pilot_radius[fold] - min(radii)) / (max(radii) - min(radii))
      for (ell in seq_along(kappa_grid)) {
        index <- smallest_minimiser(
          fold_train_risk[fold, ] + kappa_grid[ell] * penalty_shape
        )
        selected <- psi[index]
        discrepancy[fold, ell] <- max(selected - target, 0)^2 +
          log(n_validation) * max(target - selected, 0)^2
      }
    }
    mean_discrepancy <- colMeans(discrepancy)
    se_discrepancy <- apply(discrepancy, 2L, sd) / sqrt(v)
    min_index <- smallest_minimiser(mean_discrepancy)
    eligible <- which(
      mean_discrepancy <= mean_discrepancy[min_index] +
        se_discrepancy[min_index] + 1e-14
    )
    selected_kappa_index <- max(eligible)
    kappa_hat <- kappa_grid[selected_kappa_index]
    kappa_endpoint <- as.integer(
      selected_kappa_index == 1L || selected_kappa_index == length(kappa_grid)
    )
    calibration_mean_min <- mean_discrepancy[min_index]
    calibration_se_at_min <- se_discrepancy[min_index]
    full_penalty_shape <- (1 + psi) * sqrt(log(n) / n)
    pan_index <- smallest_minimiser(
      full_path$risk + kappa_hat * full_penalty_shape
    )
  }

  cv_index <- smallest_minimiser(colMeans(fold_validation_loss))
  active_count <- colSums(abs(full_path$beta) > config$active_tolerance)
  monotone_active_count <- cummax(active_count)
  if (family == "gaussian") {
    bic_criterion <- full_path$risk +
      (monotone_active_count + config$bic_epsilon * psi) * log(n) / n
    bic_label <- "BIC-like"
  } else {
    bic_criterion <- full_path$risk +
      0.5 * (monotone_active_count + config$bic_epsilon * psi) * log(n) / n
    bic_label <- "BIC-active (exploratory)"
  }
  bic_index <- smallest_minimiser(bic_criterion)

  primary <- rbind(
    if (calibration_failed) failed_method_row("PanIC-CF") else
      evaluate_path_index(full_path, pan_index, "PanIC-CF", truth, test,
                          family, kappa_hat, FALSE, config),
    evaluate_path_index(full_path, bic_index, bic_label, truth, test,
                        family, NA_real_, FALSE, config),
    evaluate_path_index(full_path, cv_index, "CV", truth, test,
                        family, NA_real_, FALSE, config)
  )

  ## Oracle candidates are retained only as explicitly labelled diagnostics.
  full_penalty_shape <- (1 + psi) * sqrt(log(n) / n)
  oracle_rows <- vector("list", length(kappa_grid))
  for (ell in seq_along(kappa_grid)) {
    index <- smallest_minimiser(
      full_path$risk + kappa_grid[ell] * full_penalty_shape
    )
    oracle_rows[[ell]] <- evaluate_path_index(
      full_path, index, "PanIC-oracle-candidate", truth, test, family,
      kappa_grid[ell], FALSE, config
    )
  }
  oracle <- do.call(rbind, oracle_rows)

  elapsed <- proc.time()[["elapsed"]] - fit_start
  common <- data.frame(
    scenario_id = scenario$scenario_id[[1L]],
    family = family,
    n = n,
    rho = rho,
    replication = replication_index,
    seed = seed,
    fold_seed = fold_seed,
    test_seed = test_seed,
    true_radius = truth$true_radius,
    signal_variance = truth$signal_variance,
    stringsAsFactors = FALSE
  )
  primary <- cbind(common[rep(1L, nrow(primary)), ], primary)
  oracle <- cbind(common[rep(1L, nrow(oracle)), ], oracle)

  diagnostic <- cbind(
    common,
    data.frame(
      calibration_failed = as.integer(calibration_failed),
      pilot_failures = sum(!pilot_ok),
      pilot_above_grid = sum(pilot_radius > max(radii), na.rm = TRUE),
      pilot_reason = paste(unique(pilot_reason[nzchar(pilot_reason)]), collapse = " | "),
      fold_warning_count = sum(fold_warning_count),
      full_warning_count = length(full_path$warnings),
      kappa_hat = kappa_hat,
      kappa_endpoint = kappa_endpoint,
      calibration_mean_min = calibration_mean_min,
      calibration_se_at_min = calibration_se_at_min,
      unregularised_radius = full_path$unregularised_radius,
      true_radius_in_grid = as.integer(truth$true_radius <= max(radii)),
      max_active_radius_error = max(
        full_path$max_active_radius_error, fold_path_radius_error
      ),
      elapsed_shared_fit_seconds = elapsed,
      radius_points = length(radii),
      radius_max = max(radii),
      stringsAsFactors = FALSE
    )
  )

  list(primary = primary, oracle = oracle, diagnostic = diagnostic)
}
