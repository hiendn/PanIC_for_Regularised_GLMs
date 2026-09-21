## Prespecified configuration for the PanIC revision simulations.
## All stochastic components use an explicit per-replication seed derived from
## master_seed; see seed_ledger.csv in the generated results.

CONFIG <- list(
  master_seed = 2146092101L,
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
  sensitivity_radius_points = c(61L, 121L, 241L),
  kappa_grid = 10^seq(-2, 2, length.out = 31L),
  folds = 5L,
  bic_epsilon = 1e-6,
  glmnet_nlambda = 300L,
  glmnet_lambda_min_ratio = 1e-6,
  glmnet_threshold = 1e-10,
  glmnet_maxit = 100000L,
  interpolation_radius_tolerance = 5e-8,
  pilot_coefficient_cap = 25,
  requested_cores = 4L,
  boundary_replications = 40000L,
  subsampling_replications = 10000L,
  runtime_replications = 20L,
  runtime_warmup = 2L
)

SCENARIOS <- data.frame(
  scenario_id = c(
    "linear_iid_n500", "linear_iid_n1000",
    "logistic_iid_n500", "logistic_iid_n1000",
    "linear_ar1_n500", "linear_ar1_n1000",
    "logistic_ar1_n500", "logistic_ar1_n1000",
    "poisson_iid_n500", "poisson_iid_n1000"
  ),
  family = c(
    "gaussian", "gaussian", "binomial", "binomial",
    "gaussian", "gaussian", "binomial", "binomial",
    "poisson", "poisson"
  ),
  n = rep(c(500L, 1000L), 5L),
  rho = c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0, 0),
  stringsAsFactors = FALSE
)

radius_grid <- function(m = CONFIG$radius_points) {
  seq(CONFIG$radius_min, CONFIG$radius_max, length.out = as.integer(m))
}

replication_seed <- function(scenario_index, replication_index) {
  as.integer(CONFIG$master_seed + 100000L * as.integer(scenario_index) +
               as.integer(replication_index))
}
