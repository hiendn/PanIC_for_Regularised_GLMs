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
} else normalizePath("Simulation_Results_Analysis.R")
production_dir <- dirname(script_path)
source(file.path(production_dir, "Simulation_Config.R"))
output_name <- arg_value("--output-dir", "results")
results_dir <- if (grepl("^/", output_name)) output_name else
  file.path(production_dir, output_name)
generated_dir <- results_dir
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

read_scenario_files <- function(suffix) {
  files <- file.path(
    results_dir, paste0(SCENARIOS$scenario_id, suffix)
  )
  files <- files[file.exists(files)]
  if (!length(files)) stop("No scenario files found for suffix ", suffix)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}
mcse <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

primary <- read_scenario_files("_primary.csv")
diagnostic <- read_scenario_files("_diagnostics.csv")
calibration <- read_scenario_files("_calibration_rows.csv")

metrics <- c(
  "fp", "fn", "fpr", "fnr", "exact", "wrong",
  "signed_attained_radius_error", "test_deviance",
  "selected_grid_radius", "attained_radius",
  "selected_lower_endpoint", "selected_upper_endpoint"
)
summary_rows <- lapply(
  split(primary, interaction(primary$scenario_id, primary$method, drop = TRUE)),
  function(dat) {
    out <- dat[1L, c("scenario_id", "family", "n", "rho", "method")]
    out$attempted_replications <- nrow(dat)
    out$failed_replications <- sum(dat$method_failed)
    out$successful_replications <- sum(!dat$method_failed)
    for (metric in metrics) {
      out[[paste0(metric, "_mean")]] <- mean(dat[[metric]], na.rm = TRUE)
      out[[paste0(metric, "_mcse")]] <- mcse(dat[[metric]])
    }
    out
  }
)
simulation_summary <- do.call(rbind, summary_rows)
simulation_summary <- simulation_summary[
  order(match(simulation_summary$scenario_id, SCENARIOS$scenario_id),
        simulation_summary$method),
]
write.csv(simulation_summary, file.path(results_dir, "simulation_summary.csv"),
          row.names = FALSE)

contrast_metrics <- c(
  "fp", "fn", "fpr", "fnr", "exact", "wrong",
  "signed_attained_radius_error", "test_deviance",
  "selected_grid_radius", "attained_radius"
)
contrast_rows <- lapply(split(primary, primary$scenario_id), function(dat) {
  lhs <- dat[dat$method == "PanIC-CF" & dat$method_failed == 0L, ]
  rhs <- dat[dat$method == "CV" & dat$method_failed == 0L, ]
  paired <- merge(
    lhs, rhs, by = c("scenario_id", "replication"),
    suffixes = c("_panic", "_cv")
  )
  if (!nrow(paired)) return(NULL)
  out <- data.frame(
    scenario_id = paired$scenario_id[1L],
    family = paired$family_panic[1L], n = paired$n_panic[1L],
    rho = paired$rho_panic[1L], comparison = "PanIC-CF minus CV",
    paired_replications = nrow(paired), stringsAsFactors = FALSE
  )
  for (metric in contrast_metrics) {
    delta <- paired[[paste0(metric, "_panic")]] -
      paired[[paste0(metric, "_cv")]]
    out[[paste0(metric, "_difference")]] <- mean(delta)
    out[[paste0(metric, "_difference_mcse")]] <- mcse(delta)
  }
  out
})
paired_contrasts <- do.call(rbind, contrast_rows)
paired_contrasts <- paired_contrasts[
  order(match(paired_contrasts$scenario_id, SCENARIOS$scenario_id)),
]
write.csv(paired_contrasts, file.path(results_dir, "paired_panic_cf_vs_cv.csv"),
          row.names = FALSE)

diagnostic_rows <- lapply(split(diagnostic, diagnostic$scenario_id), function(dat) {
  data.frame(
    scenario_id = dat$scenario_id[1L], family = dat$family[1L],
    n = dat$n[1L], rho = dat$rho[1L],
    attempted_replications = nrow(dat),
    full_path_failures = sum(dat$full_path_failed),
    calibration_failures = sum(dat$calibration_failed),
    calibration_default_uses = sum(dat$default_used),
    cv_failures = sum(dat$cv_failed),
    raw_targets_below_grid = sum(dat$raw_target_below_grid, na.rm = TRUE),
    raw_targets_above_grid = sum(dat$raw_target_above_grid, na.rm = TRUE),
    projections = sum(dat$projection_count, na.rm = TRUE),
    kappa_lower_endpoint_selections = sum(dat$kappa_lower_endpoint, na.rm = TRUE),
    kappa_upper_endpoint_selections = sum(dat$kappa_upper_endpoint, na.rm = TRUE),
    kappa_q25 = quantile(dat$kappa_hat, 0.25, na.rm = TRUE),
    kappa_median = median(dat$kappa_hat, na.rm = TRUE),
    kappa_q75 = quantile(dat$kappa_hat, 0.75, na.rm = TRUE),
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
diagnostic_summary <- diagnostic_summary[
  order(match(diagnostic_summary$scenario_id, SCENARIOS$scenario_id)),
]
write.csv(diagnostic_summary, file.path(results_dir, "diagnostic_summary.csv"),
          row.names = FALSE)

## Treat each data-set replication, rather than each of its ten dependent
## calibration rows, as the Monte Carlo unit for target summaries.
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
replication_target$projected_target_bias <-
  replication_target$projected_cross_signed_radius - replication_target$true_radius
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
    projected_target_bias = mean(dat$projected_target_bias),
    projected_target_bias_mcse = mcse(dat$projected_target_bias),
    stringsAsFactors = FALSE
  )
})
target_summary <- do.call(rbind, target_rows)
target_summary <- target_summary[
  order(match(target_summary$scenario_id, SCENARIOS$scenario_id)),
]
write.csv(target_summary, file.path(results_dir, "calibration_target_summary.csv"),
          row.names = FALSE)

scenario_label <- c(
  linear_iid_n500 = "Gaussian independent, $n=500$",
  linear_iid_n2000 = "Gaussian independent, $n=2000$",
  logistic_iid_n500 = "Logistic independent, $n=500$",
  logistic_iid_n2000 = "Logistic independent, $n=2000$",
  linear_ar1_n1000 = "Gaussian AR(1), $n=1000$",
  logistic_ar1_n1000 = "Logistic AR(1), $n=1000$",
  poisson_iid_n1000 = "Poisson independent, $n=1000$"
)
method_rank <- c("PanIC-CF", "BIC-like", "BIC-active (exploratory)", "CV")
tab <- simulation_summary
tab$scenario_rank <- match(tab$scenario_id, names(scenario_label))
tab$method_rank <- match(tab$method, method_rank)
tab <- tab[order(tab$scenario_rank, tab$method_rank), ]
cell <- function(mean, se, digits = 3L) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean, se)
}
support_lines <- c(
  "\\begin{table}[H]", "\\centering",
  paste0("\\caption{Support recovery for the revised repeated-half-split ",
         "PanIC-CF procedure and the prespecified comparators. Entries are ",
         "Monte Carlo means with Monte Carlo standard errors in parentheses ",
         "over 500 attempted replications. FPR and FNR use the ten inactive ",
         "and ten active slopes, respectively; a slope is selected when ",
         "$\\lvert\\widehat\\beta_j\\rvert>10^{-8}$.}"),
  "\\label{Table: revised primary support}",
  "\\begin{adjustbox}{width=\\textwidth}", "\\scriptsize",
  "\\begin{tabular}{llrrrrrr}", "\\hline",
  paste0("Setting & Method & FP & FN & FPR & FNR & Exact & ",
         "Total support error \\\\"),
  "\\hline"
)
last <- ""
for (i in seq_len(nrow(tab))) {
  row <- tab[i, ]
  label <- if (row$scenario_id != last) scenario_label[[row$scenario_id]] else ""
  last <- row$scenario_id
  support_lines <- c(support_lines, paste0(
    label, " & ", row$method, " & ",
    cell(row$fp_mean, row$fp_mcse, 2), " & ",
    cell(row$fn_mean, row$fn_mcse, 2), " & ",
    cell(row$fpr_mean, row$fpr_mcse, 3), " & ",
    cell(row$fnr_mean, row$fnr_mcse, 3), " & ",
    cell(row$exact_mean, row$exact_mcse, 3), " & ",
    cell(row$wrong_mean, row$wrong_mcse, 2), " \\\\"
  ))
}
support_lines <- c(support_lines, "\\hline", "\\end{tabular}",
                   "\\end{adjustbox}", "\\end{table}")
writeLines(support_lines, file.path(generated_dir, "table_primary_support.tex"))

performance_lines <- c(
  "\\begin{table}[H]", "\\centering",
  paste0("\\caption{Selected and attained constraint radii and independent-test ",
         "performance. Entries are Monte Carlo means with Monte Carlo standard ",
         "errors in parentheses over 500 attempted replications. Test deviance ",
         "is evaluated on an independent sample of size 2,000. Signed radius ",
         "error is $\\lVert\\widehat\\beta\\rVert_1-C^\\star$.}"),
  "\\label{Table: revised primary performance}",
  "\\begin{adjustbox}{width=\\textwidth}", "\\scriptsize",
  "\\begin{tabular}{llrrrr}", "\\hline",
  paste0("Setting & Method & Signed radius error & $C_{\\rm sel}$ & ",
         "$\\lVert\\widehat\\beta\\rVert_1$ & Test deviance \\\\"),
  "\\hline"
)
last <- ""
for (i in seq_len(nrow(tab))) {
  row <- tab[i, ]
  label <- if (row$scenario_id != last) scenario_label[[row$scenario_id]] else ""
  last <- row$scenario_id
  performance_lines <- c(performance_lines, paste0(
    label, " & ", row$method, " & ",
    cell(row$signed_attained_radius_error_mean,
         row$signed_attained_radius_error_mcse, 3), " & ",
    cell(row$selected_grid_radius_mean, row$selected_grid_radius_mcse, 3), " & ",
    cell(row$attained_radius_mean, row$attained_radius_mcse, 3), " & ",
    cell(row$test_deviance_mean, row$test_deviance_mcse, 3), " \\\\"
  ))
}
performance_lines <- c(performance_lines, "\\hline", "\\end{tabular}",
                       "\\end{adjustbox}", "\\end{table}")
writeLines(performance_lines,
           file.path(generated_dir, "table_primary_performance.tex"))

calibration_table <- merge(
  target_summary,
  diagnostic_summary[c(
    "scenario_id", "calibration_failures", "calibration_default_uses",
    "projections", "kappa_q25", "kappa_median", "kappa_q75"
  )],
  by = "scenario_id", all.x = TRUE, sort = FALSE
)
calibration_table <- calibration_table[
  order(match(calibration_table$scenario_id, names(scenario_label))),
]
calibration_lines <- c(
  "\\begin{table}[H]", "\\centering",
  paste0("\\caption{Calibration targets, selected multipliers, and failure ",
         "diagnostics for revised PanIC-CF. The ten raw cross-signed targets ",
         "are averaged within each data-set replication before the across-",
         "replication mean and Monte Carlo standard error in parentheses are ",
         "computed over $N=500$ replications. Raw bias is the mean target ",
         "minus $C^\\star$; multiplier entries are median [first quartile, ",
         "third quartile]. Failure, prespecified-default, and projection ",
         "columns are counts across replications or calibration rows, as ",
         "applicable.}"),
  "\\label{Table: revised calibration}",
  "\\begin{adjustbox}{width=\\textwidth}", "\\scriptsize",
  "\\begin{tabular}{lrrrlrrr}", "\\hline",
  paste0("Setting & $C^\\star$ & Mean raw target & Raw bias & ",
         "$\\widehat\\kappa$ [Q1,Q3] & Failures & Defaults & ",
         "Projections \\\\"),
  "\\hline"
)
for (i in seq_len(nrow(calibration_table))) {
  row <- calibration_table[i, ]
  calibration_lines <- c(calibration_lines, paste0(
    scenario_label[[row$scenario_id]], " & ",
    sprintf("%.3f", row$true_radius), " & ",
    cell(row$raw_target_mean, row$raw_target_mean_mcse, 3), " & ",
    cell(row$raw_target_bias, row$raw_target_bias_mcse, 3), " & ",
    sprintf("%.3f [%.3f, %.3f]", row$kappa_median,
            row$kappa_q25, row$kappa_q75), " & ",
    row$calibration_failures, " & ", row$calibration_default_uses,
    " & ", row$projections, " \\\\"
  ))
}
calibration_lines <- c(
  calibration_lines, "\\hline", "\\end{tabular}",
  "\\end{adjustbox}", "\\end{table}"
)
writeLines(calibration_lines,
           file.path(generated_dir, "table_calibration.tex"))

cat("Generated summaries and manuscript tables in ", results_dir, "\n", sep = "")
