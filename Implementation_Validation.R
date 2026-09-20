#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]))
} else {
  normalizePath("Implementation_Validation.R")
}
candidate_dir <- dirname(script_path)
.libPaths(c(file.path(candidate_dir, "Rlib"), .libPaths()))
source(file.path(candidate_dir, "Method_Lock_Verification.R"))
verify_method_lock(candidate_dir)
source(file.path(candidate_dir, "Simulation_Config.R"))
source(file.path(candidate_dir, "PanIC_CF_Functions.R"))

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, passed = as.integer(isTRUE(passed)), detail = detail,
    stringsAsFactors = FALSE
  )
}

prior_env <- new.env(parent = baseenv())
sys.source(
  file.path(candidate_dir, "reference_implementation", "config.R"),
  envir = prior_env
)
prior_config <- prior_env$CONFIG
common_fields <- intersect(names(prior_config), names(CONFIG))
changed_common_fields <- common_fields[!vapply(common_fields, function(name) {
  identical(prior_config[[name]], CONFIG[[name]])
}, logical(1))]
add_check(
  "only authorized prior-configuration fields changed",
  identical(
    sort(changed_common_fields),
    sort(c("active_tolerance", "master_seed", "n_rep"))
  ),
  paste(changed_common_fields, collapse = ",")
)
add_check(
  "confirmatory replication count is locked",
  identical(CONFIG$n_rep, 1000L),
  paste0("n_rep=", CONFIG$n_rep)
)
add_check(
  "production and dedicated-smoke master seeds are locked",
  identical(CONFIG$master_seed, 2136092001L) &&
    identical(CONFIG$grid_sensitivity_master_seed, 2138092001L) &&
    identical(CONFIG$smoke_master_seed, 2140092001L) &&
    identical(CONFIG$smoke_grid_sensitivity_master_seed, 2142092001L),
  paste0(
    "main=", CONFIG$master_seed, ";grid=",
    CONFIG$grid_sensitivity_master_seed, ";smoke_main=",
    CONFIG$smoke_master_seed, ";smoke_grid=",
    CONFIG$smoke_grid_sensitivity_master_seed
  )
)

enumerate_streams <- function(master_seed, replications,
                              scenario_indices = seq_len(nrow(SCENARIOS))) {
  unlist(lapply(scenario_indices, function(scenario_index) {
    unlist(lapply(seq_len(replications), function(replication_index) {
      base <- as.integer(
        master_seed + 100000L * scenario_index + replication_index
      )
      c(
        training = base,
        calibration = as.integer(
          base + 10000L +
            1000L * (seq_len(CONFIG$calibration_half_splits) - 1L)
        ),
        cv = as.integer(base + 30000L),
        test = as.integer(base + 60000L)
      )
    }), use.names = FALSE)
  }), use.names = FALSE)
}
main_streams <- enumerate_streams(CONFIG$master_seed, CONFIG$n_rep)
grid_streams <- enumerate_streams(
  CONFIG$grid_sensitivity_master_seed,
  CONFIG$grid_sensitivity_replications, 1L
)
smoke_main_streams <- enumerate_streams(
  CONFIG$smoke_master_seed, 1000L
)
smoke_grid_streams <- enumerate_streams(
  CONFIG$smoke_grid_sensitivity_master_seed, 500L, 1L
)
prior_streams <- unlist(lapply(CONFIG$prior_master_seeds, function(seed) {
  ## One thousand replications in all seven positions is a conservative
  ## superset of every prior main, sensitivity, and development run.
  enumerate_streams(seed, 1000L)
}), use.names = FALSE)
add_check(
  "main seed streams are internally unique",
  !anyDuplicated(main_streams),
  paste0("stream_count=", length(main_streams),
         ";range=", paste(range(main_streams), collapse = ":"))
)
add_check(
  "grid seed streams are internally unique",
  !anyDuplicated(grid_streams),
  paste0("stream_count=", length(grid_streams),
         ";range=", paste(range(grid_streams), collapse = ":"))
)
add_check(
  "all production and smoke stream families are pairwise disjoint",
  {
    stream_sets <- list(
      main = main_streams, grid = grid_streams,
      smoke_main = smoke_main_streams, smoke_grid = smoke_grid_streams
    )
    pairs <- combn(names(stream_sets), 2L, simplify = FALSE)
    all(vapply(pairs, function(pair) {
      !any(stream_sets[[pair[1L]]] %in% stream_sets[[pair[2L]]])
    }, logical(1)))
  },
  paste0(
    "counts=", paste(
      c(length(main_streams), length(grid_streams),
        length(smoke_main_streams), length(smoke_grid_streams)),
      collapse = ","
    )
  )
)
add_check(
  "production and smoke streams are disjoint from all prior ranges",
  !any(main_streams %in% prior_streams) &&
    !any(grid_streams %in% prior_streams) &&
    !any(smoke_main_streams %in% prior_streams) &&
    !any(smoke_grid_streams %in% prior_streams),
  paste0("prior_superset_stream_count=", length(prior_streams))
)
add_check(
  "all enumerated seeds remain within the R integer range",
  all(c(
    main_streams, grid_streams, smoke_main_streams, smoke_grid_streams,
    prior_streams
  ) > 0L) &&
    all(c(
      main_streams, grid_streams, smoke_main_streams, smoke_grid_streams,
      prior_streams
    ) <= .Machine$integer.max),
  paste0(
    "new_range=", paste(range(c(
      main_streams, grid_streams, smoke_main_streams, smoke_grid_streams
    )), collapse = ":")
  )
)

weight_250 <- calibration_weight(250L)
weight_large <- calibration_weight(1000000L)
add_check(
  "primary calibration weight matches locked square-root log-log formula",
  isTRUE(all.equal(
    weight_250, sqrt(log(log(250 + exp(exp(1)))))
  )),
  sprintf("omega_250=%.12f", weight_250)
)
add_check(
  "primary calibration weight has the required qualitative rates",
  weight_large > weight_250 && weight_large / 1000000 < weight_250 / 250,
  sprintf("omega_250=%.6f;omega_1e6=%.6f", weight_250, weight_large)
)

cv_loss_test <- cbind(
  rep(1.1, 5L), rep(1.005, 5L),
  c(0.98, 0.99, 1.00, 1.01, 1.02), rep(1.02, 5L)
)
cv_rule_test <- select_cv(cv_loss_test, 0:3)
add_check(
  "CV selects the smallest minimum-loss radius",
  cv_rule_test$index == 3L,
  paste0("selected_index=", cv_rule_test$index)
)

literal_support_test <- support_metrics(
  c(1, 1e-12, 0),
  list(support = c(TRUE, FALSE, FALSE), true_radius = 1)
)
add_check(
  "support uses literal fitted nonzeros without a numerical threshold",
  literal_support_test$fp == 1L && literal_support_test$fn == 0L &&
    literal_support_test$wrong == 1L,
  paste0("fp=", literal_support_test$fp, ";fn=", literal_support_test$fn)
)
add_check(
  "BIC-like active counts use literal fitted nonzeros",
  identical(
    as.integer(colSums(
      matrix(c(0, 1e-12, 0, 0, 1, -1e-15), nrow = 2L) != 0
    )),
    c(1L, 0L, 2L)
  ),
  "synthetic active counts are 1,0,2"
)

production_master_seed <- CONFIG$master_seed
CONFIG$master_seed <- CONFIG$smoke_master_seed
streams <- seed_streams(1L, 1L)
assignment <- balanced_folds(500L, 2L, streams$calibration_splits[1L])
add_check(
  "balanced half split and distinct within-replication streams",
  identical(as.integer(table(assignment)), c(250L, 250L)) &&
    !anyDuplicated(c(
      streams$training, streams$calibration_splits, streams$cv_folds,
      streams$test
    )),
  paste(names(table(assignment)), table(assignment), collapse = ";")
)

## End-to-end determinism check and a selection-only comparison with ordinary
## CV from the preceding implementation on identical data.
scenario <- SCENARIOS[1L, , drop = FALSE]
candidate_fit_1 <- run_candidate_replication(scenario, 1L, 1L)
candidate_fit_2 <- run_candidate_replication(scenario, 1L, 1L)
prior_fit <- run_revised_replication(scenario, 1L, 1L)

diagnostic_columns <- setdiff(
  names(candidate_fit_1$diagnostic), "elapsed_seconds"
)
add_check(
  "complete candidate replication is deterministic",
  isTRUE(all.equal(candidate_fit_1$primary, candidate_fit_2$primary,
                   check.attributes = TRUE)) &&
    isTRUE(all.equal(candidate_fit_1$calibration,
                     candidate_fit_2$calibration,
                     check.attributes = TRUE)) &&
    isTRUE(all.equal(
      candidate_fit_1$diagnostic[diagnostic_columns],
      candidate_fit_2$diagnostic[diagnostic_columns],
      check.attributes = TRUE
    )),
  "two executions with identical explicit streams agree"
)
expected_methods <- c(
  "PanIC-CF", "BIC-like", CONFIG$primary_cv_method
)
add_check(
  "complete replication contains exactly the three assessed methods",
  identical(candidate_fit_1$primary$method, expected_methods) &&
    all(candidate_fit_1$primary$method_failed == 0L),
  paste(candidate_fit_1$primary$method, collapse = ",")
)

add_check(
  "CV preserves the preceding ordinary minimum-loss selection",
  {
    candidate_cv <- candidate_fit_1$primary[
      candidate_fit_1$primary$method == CONFIG$primary_cv_method,
    ]
    prior_cv <- prior_fit$primary[prior_fit$primary$method == "CV", ]
    isTRUE(all.equal(
      candidate_cv$selected_grid_radius,
      prior_cv$selected_grid_radius,
      tolerance = 1e-12
    )) && isTRUE(all.equal(
      candidate_cv$test_deviance, prior_cv$test_deviance,
      tolerance = 1e-12
    ))
  },
  "selected radius and test deviance agree on identical data and folds"
)
add_check(
  "ten calibration rows retain the square-root log-log weight",
  nrow(candidate_fit_1$calibration) == 10L &&
    isTRUE(all.equal(
      candidate_fit_1$calibration$weight,
      sqrt(log(log(
        candidate_fit_1$calibration$n_validation + exp(exp(1))
      ))),
      tolerance = 1e-12, check.attributes = FALSE
    )),
  paste0("rows=", nrow(candidate_fit_1$calibration))
)
CONFIG$master_seed <- production_master_seed
add_check(
  "confirmatory decision constants are locked",
  CONFIG$primary_support_alpha == 0.05 &&
    CONFIG$prediction_noninferiority_alpha == 0.05 &&
    CONFIG$prediction_noninferiority_margin == 0.001 &&
    identical(
      CONFIG$prediction_sensitivity_margins,
      c(0.0005, 0.0025, 0.005)
    ),
  paste0(
    "alpha=", CONFIG$primary_support_alpha,
    ";primary_margin=", CONFIG$prediction_noninferiority_margin,
    ";sensitivity=",
    paste(CONFIG$prediction_sensitivity_margins, collapse = ",")
  )
)

validation <- do.call(rbind, checks)
dir.create(file.path(candidate_dir, "results"), showWarnings = FALSE)
write.csv(validation,
          file.path(candidate_dir, "results", "implementation_validation.csv"),
          row.names = FALSE)
print(validation, row.names = FALSE)
if (any(validation$passed != 1L)) {
  stop("One or more deterministic implementation checks failed", call. = FALSE)
}
cat("All confirmatory implementation checks passed.\n")
