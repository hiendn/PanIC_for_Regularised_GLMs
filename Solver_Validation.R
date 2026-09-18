#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "Solver_Validation.R"
replication_dir <- normalizePath(dirname(script_path))
source(file.path(replication_dir, "Auxiliary_Config.R"))
source(file.path(replication_dir, "Auxiliary_Simulation_Functions.R"))
results_dir <- file.path(replication_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

fit_glmnet_at_lambda <- function(x, y, family, lambda, lambda_anchor) {
  lambda_sequence <- sort(unique(c(lambda_anchor, lambda)), decreasing = TRUE)
  fit <- glmnet(
    x = x, y = y, family = family, alpha = 1,
    intercept = TRUE, standardize = FALSE,
    lambda = lambda_sequence,
    thresh = CONFIG$glmnet_threshold,
    maxit = CONFIG$glmnet_maxit
  )
  index <- which.min(abs(fit$lambda - lambda))
  list(
    beta0 = as.numeric(fit$a0[index]),
    beta = as.numeric(as.matrix(fit$beta)[, index]),
    lambda = fit$lambda[index]
  )
}

direct_radius_root <- function(x, y, family, target, tolerance = 1e-7) {
  base <- glmnet(
    x = x, y = y, family = family, alpha = 1,
    intercept = TRUE, standardize = FALSE,
    nlambda = 200L, lambda.min.ratio = 1e-7,
    thresh = CONFIG$glmnet_threshold,
    maxit = CONFIG$glmnet_maxit
  )
  beta <- as.matrix(base$beta)
  radii <- colSums(abs(beta))
  lambda <- base$lambda
  hi_index <- which(radii >= target)[1L]
  if (is.na(hi_index) || hi_index <= 1L) stop("Direct root target not bracketed")
  lo_index <- hi_index - 1L
  log_lambda_high <- log(lambda[lo_index]) # higher penalty, smaller radius
  log_lambda_low <- log(lambda[hi_index])  # lower penalty, larger radius
  anchor <- lambda[1L]
  current <- NULL
  for (iteration in seq_len(60L)) {
    log_mid <- (log_lambda_low + log_lambda_high) / 2
    current <- fit_glmnet_at_lambda(x, y, family, exp(log_mid), anchor)
    radius <- sum(abs(current$beta))
    if (abs(radius - target) <= tolerance * max(1, target)) break
    if (radius < target) {
      log_lambda_high <- log_mid
    } else {
      log_lambda_low <- log_mid
    }
  }
  current$radius <- sum(abs(current$beta))
  current$iterations <- iteration
  current
}

validation_rows <- list()
cursor <- 1L
families <- c("gaussian", "binomial", "poisson")
for (family_index in seq_along(families)) {
  family <- families[family_index]
  truth <- scenario_truth(family, 0, CONFIG)
  dat <- generate_dataset(
    500L, family, truth,
    seed = CONFIG$master_seed + 900000L + family_index,
    config = CONFIG
  )
  unregularised <- fit_unregularised(dat$x, dat$y, family, CONFIG)
  if (!unregularised$ok) stop("Validation pilot failed: ", unregularised$reason)
  targets <- round(unregularised$radius * c(0.25, 0.50, 0.75), 6)
  interpolated <- fit_radius_path(dat$x, dat$y, family, targets, CONFIG)
  if (!interpolated$ok) stop(interpolated$reason)
  for (j in seq_along(targets)) {
    direct <- direct_radius_root(dat$x, dat$y, family, targets[j])
    eta_interpolated <- drop(interpolated$beta0[j] + dat$x %*% interpolated$beta[, j])
    eta_direct <- drop(direct$beta0 + dat$x %*% direct$beta)
    validation_rows[[cursor]] <- data.frame(
      family = family,
      target_radius = targets[j],
      interpolated_radius = interpolated$attained_radius[j],
      direct_radius = direct$radius,
      max_abs_coefficient_difference = max(abs(interpolated$beta[, j] - direct$beta)),
      intercept_difference = abs(interpolated$beta0[j] - direct$beta0),
      linear_predictor_rmse = sqrt(mean((eta_interpolated - eta_direct)^2)),
      empirical_risk_difference = abs(
        mean_loss(family, dat$y, eta_interpolated) -
          mean_loss(family, dat$y, eta_direct)
      ),
      direct_root_iterations = direct$iterations,
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
  }
}
validation <- do.call(rbind, validation_rows)
write.csv(validation, file.path(results_dir, "solver_validation.csv"), row.names = FALSE)

## Separate Gamma check.  The legacy attachment used a fixed shape and did not
## implement joint shape estimation.  We therefore verify only the regression
## path numerics and do not use Gamma results in the manuscript comparison.
set.seed(CONFIG$master_seed + 910000L)
gamma_truth <- scenario_truth("gaussian", 0, CONFIG)
gamma_x <- matrix(rnorm(500L * CONFIG$d), nrow = 500L) %*% chol(gamma_truth$sigma_x)
gamma_eta <- drop(gamma_truth$beta0 + gamma_x %*% gamma_truth$beta)
gamma_shape <- 2
gamma_y <- rgamma(500L, shape = gamma_shape, scale = exp(gamma_eta) / gamma_shape)
gamma_glm <- glm.fit(
  x = cbind(1, gamma_x), y = gamma_y,
  family = Gamma(link = "log"), intercept = FALSE,
  control = glm.control(epsilon = 1e-10, maxit = 100L)
)
gamma_warnings <- character()
gamma_path <- withCallingHandlers(
  glmnet(
    x = gamma_x, y = gamma_y, family = Gamma(link = "log"), alpha = 1,
    intercept = TRUE, standardize = FALSE,
    nlambda = 300L, lambda.min.ratio = 1e-7,
    thresh = CONFIG$glmnet_threshold, maxit = CONFIG$glmnet_maxit
  ),
  warning = function(w) {
    gamma_warnings <<- c(gamma_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
gamma_last <- c(
  gamma_path$a0[length(gamma_path$a0)],
  as.numeric(as.matrix(gamma_path$beta)[, ncol(gamma_path$beta)])
)
gamma_check <- data.frame(
  fixed_shape_used_for_generation = gamma_shape,
  glm_converged = isTRUE(gamma_glm$converged),
  glm_rank = gamma_glm$rank,
  glmnet_jerr = gamma_path$jerr,
  warning_count = length(gamma_warnings),
  warning_text = paste(unique(gamma_warnings), collapse = " | "),
  smallest_lambda = min(gamma_path$lambda),
  coefficient_rmse_to_unregularised_glm = sqrt(mean((gamma_last - gamma_glm$coefficients)^2)),
  max_abs_coefficient_difference = max(abs(gamma_last - gamma_glm$coefficients)),
  stringsAsFactors = FALSE
)
write.csv(gamma_check, file.path(results_dir, "gamma_validation.csv"), row.names = FALSE)

environment_lines <- c(
  paste0("R version: ", R.version.string),
  paste0("Platform: ", R.version$platform),
  paste0("glmnet version: ", as.character(packageVersion("glmnet"))),
  paste0("Matrix version: ", as.character(packageVersion("Matrix"))),
  paste0("BLAS: ", extSoftVersion()[["BLAS"]]),
  paste0("Requested simulation workers: ", CONFIG$requested_cores),
  paste0("Active-set tolerance: ", format(CONFIG$active_tolerance, scientific = TRUE)),
  paste0("glmnet convergence threshold: ", format(CONFIG$glmnet_threshold, scientific = TRUE)),
  paste0("glmnet maximum iterations: ", CONFIG$glmnet_maxit),
  paste0("glmnet path points requested: ", CONFIG$glmnet_nlambda),
  paste0("glmnet lambda.min.ratio: ", CONFIG$glmnet_lambda_min_ratio)
)
writeLines(environment_lines, file.path(results_dir, "auxiliary_environment.txt"))
capture.output(sessionInfo(), file = file.path(results_dir, "auxiliary_sessionInfo.txt"))
