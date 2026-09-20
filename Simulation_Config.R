## Locked configuration for the confirmatory PanIC-CF simulation.
##
## The preceding implementation is vendored under reference_implementation/
## so that the scenarios, grids, numerical tolerances, and seed-stream
## construction remain identical in a fresh checkout. Only the explicitly
## listed confirmatory amendments below are permitted.

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
CONFIG$master_seed <- 2136092001L
CONFIG$n_rep <- 1000L
CONFIG$calibration_weight_power <- 0.5
CONFIG$primary_support_alpha <- 0.05
CONFIG$prediction_noninferiority_alpha <- 0.05
CONFIG$prediction_noninferiority_margin <- 0.001
CONFIG$prediction_sensitivity_margins <- c(0.0005, 0.0025, 0.005)
CONFIG$prior_master_seeds <- c(
  2026091801L, 2026091802L, 2036091802L, 2046091802L,
  2056091802L, 2066091802L, 2076091802L, 2086091802L,
  2106091802L, 2126091802L
)
CONFIG$grid_sensitivity_master_seed <- 2138092001L
CONFIG$grid_sensitivity_replications <- 500L
CONFIG$grid_sensitivity_points <- c(61L, 121L, 241L)
CONFIG$smoke_master_seed <- 2140092001L
CONFIG$smoke_grid_sensitivity_master_seed <- 2142092001L
## The assessed methods are PanIC-CF, ordinary minimum-loss five-fold CV,
## and the family-appropriate BIC-like comparator.  Support is defined by
## literal nonzero fitted coefficients, so the inherited numerical support
## threshold is set to zero for compatibility with legacy interfaces; active
## support logic uses explicit literal-nonzero comparisons.
CONFIG$primary_cv_method <- "CV"
CONFIG$active_tolerance <- 0
CONFIG$default_output_dir <- "results"

rm(candidate_config_dir, prior_production_dir)
