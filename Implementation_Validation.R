#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]))
} else normalizePath("Implementation_Validation.R")
production_dir <- dirname(script_path)
args <- commandArgs(trailingOnly = TRUE)
output_arg <- grep("^--output-dir=", args, value = TRUE)
output_name <- if (length(output_arg)) {
  sub("^--output-dir=", "", output_arg[[1L]])
} else "results"
results_dir <- if (grepl("^/", output_name)) output_name else
  file.path(production_dir, output_name)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(production_dir, "Simulation_Config.R"))
source(file.path(production_dir, "PanIC_CF_Functions.R"))

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, passed = as.integer(isTRUE(passed)), detail = detail,
    stringsAsFactors = FALSE
  )
}

streams <- seed_streams(1L, 1L)
all_streams <- c(
  streams$training, streams$calibration_splits, streams$cv_folds,
  streams$test
)
add_check(
  "seed streams are distinct",
  !anyDuplicated(all_streams),
  paste(all_streams, collapse = ",")
)

assignment <- balanced_folds(500L, 2L, streams$calibration_splits[1L])
add_check(
  "balanced half split",
  identical(as.integer(table(assignment)), c(250L, 250L)),
  paste(names(table(assignment)), table(assignment), collapse = ";")
)
first_training <- which(assignment != 1L)
first_validation <- which(assignment == 1L)
second_training <- which(assignment != 2L)
second_validation <- which(assignment == 2L)
add_check(
  "both directions exchange the halves",
  identical(first_training, second_validation) &&
    identical(first_validation, second_training) &&
    !length(intersect(first_training, first_validation)) &&
    length(union(first_training, first_validation)) == 500L,
  "direction 1 and direction 2 use complementary training and validation sets"
)

signed_test <- cross_signed_radius(c(2, -3, 0), c(4, -5, 100))
add_check(
  "cross-signed target uses sign(0)=0",
  identical(signed_test, 9),
  paste0("computed=", signed_test, "; expected=9")
)
projection_test <- project_radius(c(-1, 7, 25), 0, 20)
add_check(
  "projection is onto the closed radius interval",
  identical(as.numeric(projection_test), c(0, 7, 20)),
  paste(projection_test, collapse = ",")
)

weight_250 <- calibration_weight(250L)
weight_large <- calibration_weight(1000000L)
add_check(
  "calibration weight matches locked formula",
  isTRUE(all.equal(weight_250, log(log(250 + exp(exp(1)))))),
  sprintf("omega_250=%.12f", weight_250)
)
add_check(
  "calibration weight has the required qualitative rates",
  weight_large > weight_250 && weight_large / 1000000 < weight_250 / 250,
  sprintf("omega_250=%.6f; omega_1e6=%.6f", weight_250, weight_large)
)

one_se_test <- select_largest_kappa_one_se(
  matrix(0, nrow = 10L, ncol = 4L), c(0.01, 0.1, 1, 10)
)
add_check(
  "prespecified calibration default lies on the kappa grid",
  min(abs(CONFIG$kappa_grid - CONFIG$calibration_default_kappa)) < 1e-12 &&
    CONFIG$calibration_default_kappa == 1,
  paste0("default_kappa=", CONFIG$calibration_default_kappa)
)
add_check(
  "one-standard-error rule selects the largest eligible kappa",
  one_se_test$index == 4L && one_se_test$kappa == 10,
  paste0("selected_index=", one_se_test$index,
         "; selected_kappa=", one_se_test$kappa)
)

cv_assignment <- balanced_folds(500L, CONFIG$cv_folds, streams$cv_folds)
add_check(
  "ordinary five-fold CV allocation is balanced",
  length(table(cv_assignment)) == 5L &&
    all(as.integer(table(cv_assignment)) == 100L),
  paste(table(cv_assignment), collapse = ",")
)

## Run one complete replication twice to test the entire deterministic pipeline.
scenario <- SCENARIOS[1L, , drop = FALSE]
fit_1 <- run_revised_replication(scenario, 1L, 1L)
fit_2 <- run_revised_replication(scenario, 1L, 1L)
diagnostic_columns <- setdiff(names(fit_1$diagnostic), "elapsed_seconds")
reproducible <- isTRUE(all.equal(fit_1$primary, fit_2$primary,
                                 check.attributes = TRUE)) &&
  isTRUE(all.equal(fit_1$calibration, fit_2$calibration,
                   check.attributes = TRUE)) &&
  isTRUE(all.equal(fit_1$diagnostic[diagnostic_columns],
                   fit_2$diagnostic[diagnostic_columns],
                   check.attributes = TRUE))
add_check(
  "complete replication is deterministic",
  reproducible,
  "two executions with the same explicit seed streams agree"
)
add_check(
  "complete replication has ten calibration rows",
  nrow(fit_1$calibration) == 10L &&
    all(table(fit_1$calibration$split_repeat) == 2L) &&
    all(table(fit_1$calibration$direction) == 5L),
  paste0("rows=", nrow(fit_1$calibration))
)
add_check(
  "raw and projected targets are both retained",
  all(is.finite(fit_1$calibration$raw_cross_signed_radius)) &&
    all(is.finite(fit_1$calibration$projected_cross_signed_radius)) &&
    all(fit_1$calibration$projected_cross_signed_radius >= CONFIG$radius_min) &&
    all(fit_1$calibration$projected_cross_signed_radius <= CONFIG$radius_max),
  paste0("raw_range=", paste(range(
    fit_1$calibration$raw_cross_signed_radius
  ), collapse = ","))
)
add_check(
  "PanIC-CF and five-fold CV both complete",
  fit_1$diagnostic$calibration_failed == 0L &&
    fit_1$diagnostic$default_used == 0L &&
    fit_1$diagnostic$cv_failed == 0L &&
    all(fit_1$primary$method_failed == 0L),
  paste(fit_1$primary$method, fit_1$primary$method_failed,
        collapse = ";")
)

validation <- do.call(rbind, checks)
write.csv(validation, file.path(results_dir, "verification_checks.csv"),
          row.names = FALSE)
print(validation, row.names = FALSE)
if (any(validation$passed != 1L)) {
  stop("One or more deterministic implementation checks failed", call. = FALSE)
}
cat("All deterministic implementation checks passed.\n")
