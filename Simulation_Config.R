## Locked configuration for the confirmatory PanIC-CF simulation.
##
## The preceding implementation is vendored under reference_implementation/
## to preserve the historical implementation and configuration lineage in a
## fresh checkout. Only the explicitly listed current amendments below are
## permitted.

candidate_config_dir <- local({
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
  file.path(candidate_config_dir, "reference_implementation"),
  mustWork = TRUE
)
source(file.path(prior_production_dir, "config.R"))

## Confirmatory amendments.
CONFIG$master_seed <- 2146092101L
CONFIG$n_rep <- 1000L
CONFIG$calibration_weight_power <- 0.5
CONFIG$primary_support_alpha <- 0.05
CONFIG$prediction_noninferiority_alpha <- 0.05
CONFIG$prediction_noninferiority_margin <- 0.001
CONFIG$prediction_sensitivity_margins <- c(0.0005, 0.0025, 0.005)
CONFIG$prior_master_seeds <- c(
  2026091801L, 2026091802L, 2036091802L, 2046091802L,
  2056091802L, 2066091802L, 2076091802L, 2086091802L,
  2106091802L, 2126091802L, 2136092001L
)
CONFIG$grid_sensitivity_master_seed <- 2138092001L
CONFIG$grid_sensitivity_replications <- 500L
CONFIG$grid_sensitivity_points <- c(61L, 121L, 241L)
CONFIG$smoke_master_seed <- 2140092001L
CONFIG$smoke_grid_sensitivity_master_seed <- 2142092001L
CONFIG$runtime_master_seed <- 2096092101L
## The assessed methods are PanIC-CF, ordinary minimum-loss five-fold CV,
## and the family-appropriate BIC-like comparator.  Support is defined by
## literal nonzero fitted coefficients, so the inherited numerical support
## threshold is set to zero for compatibility with legacy interfaces; active
## support logic uses explicit literal-nonzero comparisons.
CONFIG$primary_cv_method <- "CV"
CONFIG$active_tolerance <- 0
CONFIG$default_output_dir <- "results"

## Revised primary scenario collection: each model/design pair is evaluated
## at n=500 and n=1000. The historical scenario registry remains unchanged in
## reference_implementation/config.R.
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

rm(candidate_config_dir, prior_production_dir)
