## Locked configuration for the confirmatory revised PanIC-CF simulation.
##
## This file deliberately uses a fresh master seed.  The procedure and all
## tuning grids are fixed here before the confirmatory results are inspected.

CONFIG <- list(
  master_seed = 2026091802L,
  n_rep = 500L,
  n_test = 2000L,
  d = 20L,
  active = seq.int(1L, 19L, by = 2L),
  target_linear_predictor_variance = 1,
  gaussian_sigma = 1,
  poisson_marginal_mean = 2,
  radius_min = 0,
  radius_max = 20,
  radius_points = 121L,
  kappa_grid = 10^seq(-2, 2, length.out = 31L),
  calibration_default_kappa = 1,
  calibration_half_splits = 5L,
  calibration_directions = 2L,
  cv_folds = 5L,
  active_tolerance = 1e-8,
  bic_epsilon = 1e-6,
  glmnet_nlambda = 300L,
  glmnet_lambda_min_ratio = 1e-6,
  glmnet_threshold = 1e-10,
  glmnet_maxit = 100000L,
  interpolation_radius_tolerance = 5e-8,
  pilot_coefficient_cap = 25,
  requested_cores = 4L,
  runtime_n = 1000L,
  runtime_warmup = 2L,
  runtime_replications = 20L
)

SCENARIOS <- data.frame(
  scenario_id = c(
    "linear_iid_n500", "linear_iid_n2000",
    "logistic_iid_n500", "logistic_iid_n2000",
    "linear_ar1_n1000", "logistic_ar1_n1000",
    "poisson_iid_n1000"
  ),
  family = c(
    "gaussian", "gaussian", "binomial", "binomial",
    "gaussian", "binomial", "poisson"
  ),
  n = c(500L, 2000L, 500L, 2000L, 1000L, 1000L, 1000L),
  rho = c(0, 0, 0, 0, 0.5, 0.5, 0),
  stringsAsFactors = FALSE
)

radius_grid <- function(m = CONFIG$radius_points) {
  seq(CONFIG$radius_min, CONFIG$radius_max, length.out = as.integer(m))
}

replication_seed <- function(scenario_index, replication_index) {
  as.integer(CONFIG$master_seed + 100000L * as.integer(scenario_index) +
               as.integer(replication_index))
}

seed_streams <- function(scenario_index, replication_index) {
  base <- replication_seed(scenario_index, replication_index)
  list(
    training = base,
    calibration_splits = as.integer(base + 10000L +
      1000L * (seq_len(CONFIG$calibration_half_splits) - 1L)),
    cv_folds = as.integer(base + 30000L),
    test = as.integer(base + 60000L)
  )
}
