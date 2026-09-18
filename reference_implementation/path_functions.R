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

