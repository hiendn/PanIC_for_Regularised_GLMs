#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", prefix, "="), "", hit[[1L]])
}
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]))
} else {
  normalizePath("Simulation_Results_Analysis.R")
}
candidate_dir <- dirname(script_path)
source(file.path(candidate_dir, "Method_Lock_Verification.R"))
verify_method_lock(candidate_dir)
source(file.path(candidate_dir, "Simulation_Config.R"))
output_name <- arg_value("--output-dir", CONFIG$default_output_dir)
results_dir <- if (grepl("^/", output_name)) {
  output_name
} else {
  file.path(candidate_dir, output_name)
}
generated_dir <- file.path(results_dir, "manuscript_generated")
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

read_all <- function(suffix) {
  paths <- file.path(results_dir, paste0(SCENARIOS$scenario_id, suffix))
  missing <- SCENARIOS$scenario_id[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Missing ", suffix, " files for: ", paste(missing, collapse = ", ")
    )
  }
  do.call(rbind, lapply(paths, read.csv, stringsAsFactors = FALSE))
}
mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

primary <- read_all("_primary.csv")
diagnostic <- read_all("_diagnostics.csv")
calibration <- read_all("_calibration_rows.csv")
run_configuration <- read.csv(
  file.path(results_dir, "configuration.csv"), stringsAsFactors = FALSE
)
run_n_rep <- as.integer(run_configuration$value[
  run_configuration$item == "replications_requested"
])
if (length(run_n_rep) != 1L || !is.finite(run_n_rep) || run_n_rep < 1L) {
  stop("Cannot recover the run-specific replication count")
}

metrics <- c(
  "fp", "fn", "fpr", "fnr", "exact", "wrong",
  "signed_attained_radius_error", "test_deviance",
  "selected_grid_radius", "attained_radius", "selected_lower_endpoint",
  "selected_upper_endpoint"
)
summary_rows <- lapply(
  split(primary, interaction(primary$scenario_id, primary$method, drop = TRUE)),
  function(dat) {
    out <- dat[1L, c("scenario_id", "family", "n", "rho", "method")]
    out$attempted_replications <- nrow(dat)
    out$failed_replications <- sum(dat$method_failed)
    out$successful_replications <- sum(dat$method_failed == 0L)
    for (metric in metrics) {
      out[[paste0(metric, "_mean")]] <- mean(dat[[metric]], na.rm = TRUE)
      out[[paste0(metric, "_mcse")]] <- mcse(dat[[metric]])
    }
    out
  }
)
simulation_summary <- do.call(rbind, summary_rows)
method_order <- c(
  "PanIC-CF", CONFIG$original_sensitivity_method, "BIC-like",
  "BIC-active (exploratory)", CONFIG$primary_cv_method,
  CONFIG$secondary_cv_method
)
simulation_summary <- simulation_summary[order(
  match(simulation_summary$scenario_id, SCENARIOS$scenario_id),
  match(simulation_summary$method, method_order)
), ]
write.csv(simulation_summary, file.path(results_dir, "simulation_summary.csv"),
          row.names = FALSE)

contrast_metrics <- c(
  "fp", "fn", "fpr", "fnr", "exact", "wrong",
  "signed_attained_radius_error", "test_deviance",
  "selected_grid_radius", "attained_radius"
)
method_pairs <- list(
  c("PanIC-CF", CONFIG$primary_cv_method),
  c(CONFIG$original_sensitivity_method, CONFIG$primary_cv_method),
  c("PanIC-CF", CONFIG$original_sensitivity_method),
  c("PanIC-CF", CONFIG$secondary_cv_method),
  c(CONFIG$original_sensitivity_method, CONFIG$secondary_cv_method),
  c(CONFIG$secondary_cv_method, CONFIG$primary_cv_method)
)
paired_method_contrast <- function(dat, lhs_method, rhs_method) {
  lhs <- dat[dat$method == lhs_method & dat$method_failed == 0L, ]
  rhs <- dat[dat$method == rhs_method & dat$method_failed == 0L, ]
  paired <- merge(
    lhs, rhs, by = c("scenario_id", "replication"),
    suffixes = c("_lhs", "_rhs")
  )
  if (!nrow(paired)) return(NULL)
  out <- data.frame(
    scenario_id = paired$scenario_id[1L], family = paired$family_lhs[1L],
    n = paired$n_lhs[1L], rho = paired$rho_lhs[1L],
    lhs_method = lhs_method, rhs_method = rhs_method,
    comparison = paste(lhs_method, "minus", rhs_method),
    paired_replications = nrow(paired), stringsAsFactors = FALSE
  )
  for (metric in contrast_metrics) {
    delta <- paired[[paste0(metric, "_lhs")]] -
      paired[[paste0(metric, "_rhs")]]
    out[[paste0(metric, "_difference")]] <- mean(delta)
    out[[paste0(metric, "_difference_mcse")]] <- mcse(delta)
  }
  relative_deviance <-
    (paired$test_deviance_lhs - paired$test_deviance_rhs) /
    paired$test_deviance_rhs
  out$relative_test_deviance_difference <- mean(relative_deviance)
  out$relative_test_deviance_difference_mcse <- mcse(relative_deviance)
  out
}
paired_rows <- list()
row_index <- 1L
for (scenario_id in SCENARIOS$scenario_id) {
  dat <- primary[primary$scenario_id == scenario_id, , drop = FALSE]
  for (pair in method_pairs) {
    paired_rows[[row_index]] <- paired_method_contrast(dat, pair[1L], pair[2L])
    row_index <- row_index + 1L
  }
}
paired_contrasts <- do.call(rbind, paired_rows)
paired_contrasts <- paired_contrasts[order(
  match(paired_contrasts$scenario_id, SCENARIOS$scenario_id),
  match(
    paired_contrasts$comparison,
    vapply(method_pairs, function(pair) paste(pair[1L], "minus", pair[2L]),
           character(1))
  )
), ]
write.csv(paired_contrasts,
          file.path(results_dir, "paired_method_contrasts.csv"),
          row.names = FALSE)
write.csv(
  paired_contrasts[
    paired_contrasts$lhs_method == "PanIC-CF" &
      paired_contrasts$rhs_method == CONFIG$primary_cv_method,
  ],
  file.path(results_dir, "paired_panic_cf_vs_cv_min.csv"), row.names = FALSE
)

z_support <- qnorm(1 - CONFIG$primary_support_alpha)
z_prediction <- qnorm(1 - CONFIG$prediction_noninferiority_alpha)
primary_scenario <- paired_contrasts[
  paired_contrasts$lhs_method == "PanIC-CF" &
    paired_contrasts$rhs_method == CONFIG$primary_cv_method,
  c(
    "scenario_id", "family", "n", "rho", "paired_replications",
    "wrong_difference", "wrong_difference_mcse",
    "relative_test_deviance_difference",
    "relative_test_deviance_difference_mcse"
  )
]
primary_scenario <- primary_scenario[match(
  SCENARIOS$scenario_id, primary_scenario$scenario_id
), ]
primary_scenario$support_upper_one_sided_95 <-
  primary_scenario$wrong_difference +
  z_support * primary_scenario$wrong_difference_mcse
primary_scenario$relative_deviance_upper_one_sided_95 <-
  primary_scenario$relative_test_deviance_difference +
  z_prediction *
    primary_scenario$relative_test_deviance_difference_mcse
write.csv(primary_scenario,
          file.path(results_dir, "scenario_primary_estimands.csv"),
          row.names = FALSE)

complete_pairing <- all(
  primary_scenario$paired_replications == run_n_rep
)
support_estimate <- mean(primary_scenario$wrong_difference)
support_mcse <- sqrt(sum(primary_scenario$wrong_difference_mcse^2)) /
  nrow(SCENARIOS)
support_upper <- support_estimate + z_support * support_mcse
prediction_estimate <- mean(
  primary_scenario$relative_test_deviance_difference
)
prediction_mcse <- sqrt(sum(
  primary_scenario$relative_test_deviance_difference_mcse^2
)) / nrow(SCENARIOS)
prediction_upper <- prediction_estimate + z_prediction * prediction_mcse
support_pass <- isTRUE(is.finite(support_upper) && support_upper < 0)
prediction_pass <- isTRUE(
  is.finite(prediction_upper) &&
    prediction_upper < CONFIG$prediction_noninferiority_margin
)
confirmatory_decision <- data.frame(
  scenarios = nrow(SCENARIOS),
  replications_per_scenario = run_n_rep,
  complete_pairing = as.integer(complete_pairing),
  support_estimand =
    "equal-weight scenario mean paired total-support-error difference",
  support_estimate = support_estimate, support_mcse = support_mcse,
  support_upper_one_sided_95 = support_upper,
  support_superiority_threshold = 0,
  support_superiority_pass = as.integer(support_pass),
  prediction_estimand =
    "equal-weight scenario mean paired relative test-deviance difference",
  prediction_estimate = prediction_estimate,
  prediction_mcse = prediction_mcse,
  prediction_upper_one_sided_95 = prediction_upper,
  prediction_noninferiority_margin =
    CONFIG$prediction_noninferiority_margin,
  prediction_noninferiority_pass = as.integer(prediction_pass),
  joint_claim_pass = as.integer(
    complete_pairing && support_pass && prediction_pass
  ),
  stringsAsFactors = FALSE
)
write.csv(confirmatory_decision,
          file.path(results_dir, "confirmatory_decision.csv"),
          row.names = FALSE)

sensitivity_margins <- data.frame(
  margin = c(
    CONFIG$prediction_noninferiority_margin,
    CONFIG$prediction_sensitivity_margins
  ),
  role = c(
    "primary",
    rep("sensitivity", length(CONFIG$prediction_sensitivity_margins))
  ),
  upper_one_sided_95 = prediction_upper,
  noninferiority_pass = as.integer(
    prediction_upper < c(
      CONFIG$prediction_noninferiority_margin,
      CONFIG$prediction_sensitivity_margins
    )
  ),
  stringsAsFactors = FALSE
)
write.csv(sensitivity_margins,
          file.path(results_dir, "prediction_margin_sensitivity.csv"),
          row.names = FALSE)

diagnostic_rows <- lapply(split(diagnostic, diagnostic$scenario_id), function(dat) {
  data.frame(
    scenario_id = dat$scenario_id[1L], family = dat$family[1L],
    n = dat$n[1L], rho = dat$rho[1L], attempted_replications = nrow(dat),
    full_path_failures = sum(dat$full_path_failed),
    calibration_failures = sum(dat$calibration_failed),
    calibration_default_uses = sum(dat$default_used),
    cv_failures = sum(dat$cv_failed),
    projections = sum(dat$projection_count, na.rm = TRUE),
    kappa_lower_endpoint_selections =
      sum(dat$kappa_lower_endpoint, na.rm = TRUE),
    kappa_upper_endpoint_selections =
      sum(dat$kappa_upper_endpoint, na.rm = TRUE),
    original_kappa_lower_endpoint_selections =
      sum(dat$original_kappa_lower_endpoint, na.rm = TRUE),
    original_kappa_upper_endpoint_selections =
      sum(dat$original_kappa_upper_endpoint, na.rm = TRUE),
    kappa_q25 = quantile(dat$kappa_hat, 0.25, na.rm = TRUE),
    kappa_median = median(dat$kappa_hat, na.rm = TRUE),
    kappa_q75 = quantile(dat$kappa_hat, 0.75, na.rm = TRUE),
    original_kappa_median = median(dat$original_kappa_hat, na.rm = TRUE),
    solver_warnings = sum(
      dat$full_warning_count + dat$calibration_warning_count +
        dat$cv_warning_count, na.rm = TRUE
    ),
    maximum_radius_interpolation_error = max(
      dat$maximum_active_radius_error, na.rm = TRUE
    ),
    mean_replication_seconds = mean(dat$elapsed_seconds, na.rm = TRUE),
    mcse_replication_seconds = mcse(dat$elapsed_seconds),
    stringsAsFactors = FALSE
  )
})
diagnostic_summary <- do.call(rbind, diagnostic_rows)
diagnostic_summary <- diagnostic_summary[match(
  SCENARIOS$scenario_id, diagnostic_summary$scenario_id
), ]
write.csv(diagnostic_summary, file.path(results_dir, "diagnostic_summary.csv"),
          row.names = FALSE)

calibration_success <- calibration[
  is.finite(calibration$raw_cross_signed_radius),
]
replication_target <- aggregate(
  cbind(raw_cross_signed_radius, projected_cross_signed_radius) ~
    scenario_id + family + n + rho + replication + true_radius,
  calibration_success, mean
)
replication_target$raw_target_bias <-
  replication_target$raw_cross_signed_radius - replication_target$true_radius
target_rows <- lapply(split(replication_target, replication_target$scenario_id),
                      function(dat) {
  data.frame(
    scenario_id = dat$scenario_id[1L], family = dat$family[1L],
    n = dat$n[1L], rho = dat$rho[1L], true_radius = dat$true_radius[1L],
    complete_replications = nrow(dat),
    raw_target_mean = mean(dat$raw_cross_signed_radius),
    raw_target_mean_mcse = mcse(dat$raw_cross_signed_radius),
    raw_target_bias = mean(dat$raw_target_bias),
    raw_target_bias_mcse = mcse(dat$raw_target_bias),
    projected_target_mean = mean(dat$projected_cross_signed_radius),
    projected_target_mean_mcse = mcse(dat$projected_cross_signed_radius),
    stringsAsFactors = FALSE
  )
})
target_summary <- do.call(rbind, target_rows)
target_summary <- target_summary[match(
  SCENARIOS$scenario_id, target_summary$scenario_id
), ]
write.csv(target_summary,
          file.path(results_dir, "calibration_target_summary.csv"),
          row.names = FALSE)

scenario_labels <- c(
  linear_iid_n500 = "Gaussian independent, $n=500$",
  linear_iid_n2000 = "Gaussian independent, $n=2000$",
  logistic_iid_n500 = "Logistic independent, $n=500$",
  logistic_iid_n2000 = "Logistic independent, $n=2000$",
  linear_ar1_n1000 = "Gaussian AR(1), $n=1000$",
  logistic_ar1_n1000 = "Logistic AR(1), $n=1000$",
  poisson_iid_n1000 = "Poisson independent, $n=1000$"
)
cell <- function(estimate, se, digits = 3L) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), estimate, se)
}
reader_method <- function(method) {
  ifelse(method == "BIC-active (exploratory)", "BIC-like", method)
}
support_lines <- c(
  "\\begin{table}[H]", "\\centering",
  paste0(
    "\\caption{Support recovery for PanIC-CF and the prespecified ",
    "comparators. ",
    "Entries are Monte Carlo means with Monte Carlo standard errors in ",
    "parentheses over ", run_n_rep, " attempted replications. ",
    "PanIC-CF-original is a same-path sensitivity analysis; CV-1SE is a ",
    "secondary trade-off comparator. The logistic and Poisson BIC-like ",
    "rows are exploratory active-count analogues outside the scope of the ",
    "Gaussian BIC-like proposition.}"
  ),
  "\\label{Table: second confirmation support}",
  "\\begin{adjustbox}{width=\\textwidth}", "\\scriptsize",
  "\\begin{tabular}{llrrrr}", "\\hline",
  "Setting & Method & FP & FN & Exact & Total support error \\\\",
  "\\hline"
)
last <- ""
for (i in seq_len(nrow(simulation_summary))) {
  row <- simulation_summary[i, ]
  label <- if (row$scenario_id != last) {
    scenario_labels[[row$scenario_id]]
  } else {
    ""
  }
  last <- row$scenario_id
  support_lines <- c(support_lines, paste0(
    label, " & ", reader_method(row$method), " & ",
    cell(row$fp_mean, row$fp_mcse, 2), " & ",
    cell(row$fn_mean, row$fn_mcse, 2), " & ",
    cell(row$exact_mean, row$exact_mcse, 3), " & ",
    cell(row$wrong_mean, row$wrong_mcse, 2), " \\\\"
  ))
}
support_lines <- c(
  support_lines, "\\hline", "\\end{tabular}",
  "\\end{adjustbox}", "\\end{table}"
)
writeLines(support_lines,
           file.path(generated_dir, "table_primary_support.tex"))

decision_lines <- c(
  "\\begin{table}[H]", "\\centering",
  "\\caption{Prespecified joint confirmatory decision.}",
  "\\label{Table: second confirmation decision}",
  "\\begin{tabular}{lrrrr}", "\\hline",
  "Endpoint & Estimate & MCSE & One-sided 95\\% upper bound & Criterion \\\\",
  "\\hline",
  sprintf(
    "Total support error & %.5f & %.5f & %.5f & $<0$ \\\\",
    support_estimate, support_mcse, support_upper
  ),
  sprintf(
    "Relative test deviance & %.6f & %.6f & %.6f & $<%.3f$ \\\\",
    prediction_estimate, prediction_mcse, prediction_upper,
    CONFIG$prediction_noninferiority_margin
  ),
  "\\hline", "\\end{tabular}", "\\end{table}"
)
writeLines(decision_lines,
           file.path(generated_dir, "table_confirmatory_decision.tex"))

cat("Confirmatory summaries and manuscript tables written to ", results_dir,
    ".\n", sep = "")
