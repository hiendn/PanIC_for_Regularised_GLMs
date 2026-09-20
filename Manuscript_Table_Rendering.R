## Deterministic renderers for the seven manuscript table fragments.
##
## Each render_* function is pure: it accepts already-read result objects and
## returns the complete LaTeX fragment as a character vector.  The thin I/O
## layer at the end reads the committed summary CSVs and writes UTF-8 files
## with LF line endings and one final newline.

PANIC_RENDERED_TABLE_FILES <- c(
  "table_primary_support.tex",
  "table_primary_performance.tex",
  "table_confirmatory_decision.tex",
  "table_calibration.tex",
  "table_grid_sensitivity.tex",
  "table_boundary.tex",
  "table_subsampling.tex"
)

PANIC_SCENARIO_LABELS <- c(
  linear_iid_n500 = "Gaussian independent, $n=500$",
  linear_iid_n2000 = "Gaussian independent, $n=2000$",
  logistic_iid_n500 = "Logistic independent, $n=500$",
  logistic_iid_n2000 = "Logistic independent, $n=2000$",
  linear_ar1_n1000 = "Gaussian AR(1), $n=1000$",
  logistic_ar1_n1000 = "Logistic AR(1), $n=1000$",
  poisson_iid_n1000 = "Poisson independent, $n=1000$"
)

PANIC_METHOD_ORDER <- c(
  "PanIC-CF", "PanIC-CF-original", "BIC-like", "CV-min", "CV-1SE"
)

panic_require_columns <- function(data, required, object_name) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

panic_reader_method <- function(method) {
  ifelse(method == "BIC-active (exploratory)", "BIC-like", method)
}

panic_cell <- function(estimate, mcse, digits) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), estimate, mcse)
}

panic_plain_integer <- function(value) {
  if (length(value) != 1L || !is.finite(value) || value != round(value)) {
    stop("Expected one finite integer-valued quantity", call. = FALSE)
  }
  formatC(as.integer(value), format = "d", big.mark = ",")
}

panic_latex_integer <- function(value) {
  gsub(",", "{,}", panic_plain_integer(value), fixed = TRUE)
}

panic_single_value <- function(value, description) {
  unique_value <- unique(value)
  if (length(unique_value) != 1L || is.na(unique_value)) {
    stop(description, " must have one non-missing value", call. = FALSE)
  }
  unique_value[[1L]]
}

panic_configuration_value <- function(configuration, item) {
  panic_require_columns(configuration, c("item", "value"), "configuration")
  hit <- configuration$value[configuration$item == item]
  if (length(hit) != 1L || is.na(hit)) {
    stop("Configuration item must occur exactly once: ", item,
         call. = FALSE)
  }
  hit[[1L]]
}

panic_order_primary_rows <- function(summary) {
  panic_require_columns(summary, c("scenario_id", "method"),
                        "simulation_summary")
  display_method <- panic_reader_method(summary$method)
  scenario_rank <- match(summary$scenario_id, names(PANIC_SCENARIO_LABELS))
  method_rank <- match(display_method, PANIC_METHOD_ORDER)
  if (anyNA(scenario_rank) || anyNA(method_rank)) {
    stop("Unexpected scenario or method in simulation_summary", call. = FALSE)
  }
  keys <- paste(summary$scenario_id, display_method, sep = "\r")
  expected_keys <- as.vector(outer(
    names(PANIC_SCENARIO_LABELS), PANIC_METHOD_ORDER, paste, sep = "\r"
  ))
  if (nrow(summary) != length(expected_keys) || anyDuplicated(keys) ||
      !setequal(keys, expected_keys)) {
    stop("simulation_summary does not contain the expected 7-by-5 rows",
         call. = FALSE)
  }
  ordered <- summary[order(scenario_rank, method_rank), , drop = FALSE]
  ordered$reader_method <- panic_reader_method(ordered$method)
  rownames(ordered) <- NULL
  ordered
}

render_primary_support_table <- function(simulation_summary, configuration) {
  required <- c(
    "scenario_id", "method", "attempted_replications",
    "fp_mean", "fp_mcse", "fn_mean", "fn_mcse", "fpr_mean", "fpr_mcse",
    "fnr_mean", "fnr_mcse", "exact_mean", "exact_mcse", "wrong_mean",
    "wrong_mcse"
  )
  panic_require_columns(simulation_summary, required, "simulation_summary")
  rows <- panic_order_primary_rows(simulation_summary)
  attempted <- panic_single_value(
    rows$attempted_replications, "attempted replications"
  )
  dimension <- as.numeric(panic_configuration_value(configuration, "dimension"))
  active <- as.numeric(panic_configuration_value(configuration, "active_slopes"))
  if (!identical(dimension - active, 10) || !identical(active, 10)) {
    stop("The current support-table caption requires ten active and ten inactive slopes",
         call. = FALSE)
  }

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0(
      "\\caption{Support recovery for PanIC-CF and the prespecified ",
      "comparators. Entries are Monte Carlo means with Monte Carlo standard ",
      "errors in parentheses over ", panic_plain_integer(attempted),
      " attempted replications. PanIC-CF-original is a same-path sensitivity ",
      "analysis, and CV-1SE is a secondary trade-off comparator. FPR and FNR ",
      "use the ten inactive and ten active slopes, respectively; a slope is ",
      "selected when $\\lvert\\widehat\\beta_j\\rvert>10^{-8}$. The active-count ",
      "comparator is labelled BIC-like in every family for compactness, but ",
      "its logistic and Poisson rows are exploratory analogues outside the ",
      "scope of Proposition \\ref{Prop: BIC-like consistency}.}"
    ),
    "\\label{Table: revised primary support}",
    "\\begin{adjustbox}{width=\\textwidth}",
    "\\scriptsize",
    "\\begin{tabular}{llrrrrrr}",
    "\\hline",
    paste0(
      "Setting & Method & FP & FN & FPR & FNR & Exact & Total support error ",
      "\\\\"
    ),
    "\\hline"
  )

  previous_scenario <- ""
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, , drop = FALSE]
    setting <- if (row$scenario_id != previous_scenario) {
      PANIC_SCENARIO_LABELS[[row$scenario_id]]
    } else {
      ""
    }
    previous_scenario <- row$scenario_id
    lines <- c(lines, paste0(
      setting, " & ", row$reader_method, " & ",
      panic_cell(row$fp_mean, row$fp_mcse, 2L), " & ",
      panic_cell(row$fn_mean, row$fn_mcse, 2L), " & ",
      panic_cell(row$fpr_mean, row$fpr_mcse, 3L), " & ",
      panic_cell(row$fnr_mean, row$fnr_mcse, 3L), " & ",
      panic_cell(row$exact_mean, row$exact_mcse, 3L), " & ",
      panic_cell(row$wrong_mean, row$wrong_mcse, 2L), " \\\\"
    ))
  }
  c(lines, "\\hline", "\\end{tabular}", "\\end{adjustbox}",
    "\\end{table}")
}

render_primary_performance_table <- function(simulation_summary,
                                             configuration) {
  required <- c(
    "scenario_id", "method", "attempted_replications",
    "signed_attained_radius_error_mean",
    "signed_attained_radius_error_mcse", "selected_grid_radius_mean",
    "selected_grid_radius_mcse", "attained_radius_mean",
    "attained_radius_mcse", "test_deviance_mean", "test_deviance_mcse"
  )
  panic_require_columns(simulation_summary, required, "simulation_summary")
  rows <- panic_order_primary_rows(simulation_summary)
  attempted <- panic_single_value(
    rows$attempted_replications, "attempted replications"
  )
  test_size <- as.numeric(panic_configuration_value(
    configuration, "independent_test_size"
  ))

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0(
      "\\caption{Selected and attained constraint radii and independent-test ",
      "performance. Entries are Monte Carlo means with Monte Carlo standard ",
      "errors in parentheses over ", panic_plain_integer(attempted),
      " attempted replications. Test deviance is evaluated on an independent ",
      "sample of size ", panic_plain_integer(test_size), ". Signed radius error ",
      "is $\\lVert\\widehat\\beta\\rVert_1-C^\\star$. PanIC-CF-original is a ",
      "same-path sensitivity analysis, and CV-1SE is a secondary trade-off ",
      "comparator. The active-count comparator is labelled BIC-like in every ",
      "family for compactness, but its logistic and Poisson rows are exploratory ",
      "analogues outside the scope of Proposition ",
      "\\ref{Prop: BIC-like consistency}.}"
    ),
    "\\label{Table: revised primary performance}",
    "\\begin{adjustbox}{width=\\textwidth}",
    "\\scriptsize",
    "\\begin{tabular}{llrrrr}",
    "\\hline",
    paste0(
      "Setting & Method & Signed radius error & $C_{\\rm sel}$ & ",
      "$\\lVert\\widehat\\beta\\rVert_1$ & Test deviance \\\\"
    ),
    "\\hline"
  )

  previous_scenario <- ""
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, , drop = FALSE]
    setting <- if (row$scenario_id != previous_scenario) {
      PANIC_SCENARIO_LABELS[[row$scenario_id]]
    } else {
      ""
    }
    previous_scenario <- row$scenario_id
    lines <- c(lines, paste0(
      setting, " & ", row$reader_method, " & ",
      panic_cell(
        row$signed_attained_radius_error_mean,
        row$signed_attained_radius_error_mcse, 3L
      ), " & ",
      panic_cell(
        row$selected_grid_radius_mean, row$selected_grid_radius_mcse, 3L
      ), " & ",
      panic_cell(row$attained_radius_mean, row$attained_radius_mcse, 3L),
      " & ",
      panic_cell(row$test_deviance_mean, row$test_deviance_mcse, 4L),
      " \\\\"
    ))
  }
  c(lines, "\\hline", "\\end{tabular}", "\\end{adjustbox}",
    "\\end{table}")
}

render_confirmatory_decision_table <- function(confirmatory_decision) {
  required <- c(
    "support_estimate", "support_mcse", "support_upper_one_sided_95",
    "support_superiority_threshold", "prediction_estimate",
    "prediction_mcse", "prediction_upper_one_sided_95",
    "prediction_noninferiority_margin"
  )
  panic_require_columns(confirmatory_decision, required,
                        "confirmatory_decision")
  if (nrow(confirmatory_decision) != 1L) {
    stop("confirmatory_decision must contain exactly one row", call. = FALSE)
  }
  row <- confirmatory_decision[1L, , drop = FALSE]
  if (row$support_superiority_threshold != 0 ||
      row$prediction_noninferiority_margin != 0.001) {
    stop("The current decision-table criteria have changed", call. = FALSE)
  }

  c(
    "\\begin{table}[H]",
    "\\centering",
    paste0(
      "\\caption{Prespecified joint confirmatory comparison of PanIC-CF with ",
      "CV-min. The two endpoints are equal-weight means over the seven ",
      "settings; the relative-deviance endpoint is dimensionless. The joint ",
      "conclusion requires both one-sided criteria to hold.}"
    ),
    "\\label{Table: confirmatory decision}",
    "\\begin{adjustbox}{width=\\textwidth}",
    "\\begin{tabular}{lrrrr}",
    "\\hline",
    paste0(
      "Endpoint & Estimate & MCSE & One-sided 95\\% upper bound & Criterion ",
      "\\\\"
    ),
    "\\hline",
    sprintf(
      paste0(
        "Total support error, $\\Delta_{\\rm S}$ & %.5f & %.5f & %.5f & ",
        "$<0$ \\\\"
      ),
      row$support_estimate, row$support_mcse,
      row$support_upper_one_sided_95
    ),
    sprintf(
      paste0(
        "Relative test deviance, $\\Delta_{\\rm D}$ & %.6f & %.6f & %.6f & ",
        "$<%.3f$ \\\\"
      ),
      row$prediction_estimate, row$prediction_mcse,
      row$prediction_upper_one_sided_95,
      row$prediction_noninferiority_margin
    ),
    "\\hline",
    "\\end{tabular}",
    "\\end{adjustbox}",
    "\\end{table}"
  )
}

render_calibration_table <- function(calibration_target_summary,
                                     diagnostic_summary) {
  target_required <- c(
    "scenario_id", "complete_replications", "true_radius",
    "raw_target_mean", "raw_target_mean_mcse", "raw_target_bias",
    "raw_target_bias_mcse"
  )
  diagnostic_required <- c(
    "scenario_id", "attempted_replications", "calibration_failures",
    "calibration_default_uses", "projections", "kappa_q25",
    "kappa_median", "kappa_q75"
  )
  panic_require_columns(calibration_target_summary, target_required,
                        "calibration_target_summary")
  panic_require_columns(diagnostic_summary, diagnostic_required,
                        "diagnostic_summary")
  if (anyDuplicated(calibration_target_summary$scenario_id) ||
      anyDuplicated(diagnostic_summary$scenario_id)) {
    stop("Calibration summaries must have unique scenario rows", call. = FALSE)
  }
  rows <- merge(
    calibration_target_summary,
    diagnostic_summary[, diagnostic_required, drop = FALSE],
    by = "scenario_id", sort = FALSE
  )
  rank <- match(rows$scenario_id, names(PANIC_SCENARIO_LABELS))
  if (nrow(rows) != length(PANIC_SCENARIO_LABELS) || anyNA(rank) ||
      !setequal(rows$scenario_id, names(PANIC_SCENARIO_LABELS))) {
    stop("Calibration summaries do not contain the expected seven scenarios",
         call. = FALSE)
  }
  rows <- rows[order(rank), , drop = FALSE]
  replications <- panic_single_value(
    rows$attempted_replications, "calibration attempted replications"
  )
  if (any(rows$complete_replications != replications)) {
    stop("Calibration target summary is not complete", call. = FALSE)
  }

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0(
      "\\caption{Calibration targets, selected multipliers, and failure ",
      "diagnostics for PanIC-CF. The ten raw cross-signed targets are averaged ",
      "within each data-set replication before the across-replication mean and ",
      "Monte Carlo standard error in parentheses are computed over $N=",
      panic_latex_integer(replications), "$ replications. Raw bias is the mean ",
      "target minus $C^\\star$; multiplier entries are median [first quartile, ",
      "third quartile]. Calibration-failure and prespecified-default columns ",
      "count replications, whereas projections count calibration rows.}"
    ),
    "\\label{Table: revised calibration}",
    "\\begin{adjustbox}{width=\\textwidth}",
    "\\scriptsize",
    "\\begin{tabular}{lrrrlrrr}",
    "\\hline",
    paste0(
      "Setting & $C^\\star$ & Mean raw target & Raw bias & ",
      "$\\widehat\\kappa$ [Q1,Q3] & Cal. failures & Defaults & Projections ",
      "\\\\"
    ),
    "\\hline"
  )
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, , drop = FALSE]
    lines <- c(lines, paste0(
      PANIC_SCENARIO_LABELS[[row$scenario_id]], " & ",
      sprintf("%.3f", row$true_radius), " & ",
      panic_cell(row$raw_target_mean, row$raw_target_mean_mcse, 3L), " & ",
      panic_cell(row$raw_target_bias, row$raw_target_bias_mcse, 3L), " & ",
      sprintf(
        "%.3f [%.3f, %.3f]", row$kappa_median, row$kappa_q25,
        row$kappa_q75
      ), " & ",
      sprintf("%d", as.integer(row$calibration_failures)), " & ",
      sprintf("%d", as.integer(row$calibration_default_uses)), " & ",
      sprintf("%d", as.integer(row$projections)), " \\\\"
    ))
  }
  c(lines, "\\hline", "\\end{tabular}", "\\end{adjustbox}",
    "\\end{table}")
}

render_grid_sensitivity_table <- function(grid_sensitivity_summary) {
  required <- c(
    "radius_points", "method", "grid_spacing", "attempted_replications",
    "wrong_mean", "wrong_mcse", "exact_mean", "exact_mcse",
    "test_deviance_mean", "test_deviance_mcse",
    "selected_grid_radius_mean", "selected_grid_radius_mcse",
    "mean_shared_runtime", "mcse_shared_runtime"
  )
  panic_require_columns(grid_sensitivity_summary, required,
                        "grid_sensitivity_summary")
  rows <- grid_sensitivity_summary
  rows$reader_method <- panic_reader_method(rows$method)
  radius_values <- sort(unique(rows$radius_points))
  radius_rank <- match(rows$radius_points, radius_values)
  method_rank <- match(rows$reader_method, PANIC_METHOD_ORDER)
  keys <- paste(rows$radius_points, rows$reader_method, sep = "\r")
  expected_keys <- as.vector(outer(
    radius_values, PANIC_METHOD_ORDER, paste, sep = "\r"
  ))
  if (anyNA(method_rank) || nrow(rows) != length(expected_keys) ||
      anyDuplicated(keys) || !setequal(keys, expected_keys)) {
    stop("grid_sensitivity_summary is not a complete grid-by-method table",
         call. = FALSE)
  }
  rows <- rows[order(radius_rank, method_rank), , drop = FALSE]
  replications <- panic_single_value(
    rows$attempted_replications, "grid-sensitivity attempted replications"
  )

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0(
      "\\caption{Radius-grid sensitivity in the independent Gaussian setting ",
      "with $n=1000$ and common random numbers across $N=",
      panic_plain_integer(replications), "$ replications. Every grid spans ",
      "$[0,20]$. Entries are means with Monte Carlo standard errors in ",
      "parentheses. PanIC-CF-original and PanIC-CF share fitted paths, as do ",
      "CV-min and CV-1SE. Shared runtime is the end-to-end time for the ",
      "full-sample path, ten half-sample calibration directions, and five ",
      "ordinary CV folds.}"
    ),
    "\\label{Table: grid sensitivity}",
    "\\begin{adjustbox}{width=\\textwidth}",
    "\\scriptsize",
    "\\begin{tabular}{rrlrrrrr}",
    "\\hline",
    paste0(
      "$m$ & Spacing & Method & Total support error & Exact & Test deviance ",
      "& Selected radius & Shared runtime (s) \\\\"
    ),
    "\\hline"
  )
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, , drop = FALSE]
    lines <- c(lines, paste0(
      sprintf("%d", as.integer(row$radius_points)), " & ",
      sprintf("%.3f", row$grid_spacing), " & ", row$reader_method, " & ",
      panic_cell(row$wrong_mean, row$wrong_mcse, 2L), " & ",
      panic_cell(row$exact_mean, row$exact_mcse, 3L), " & ",
      panic_cell(row$test_deviance_mean, row$test_deviance_mcse, 3L),
      " & ",
      panic_cell(
        row$selected_grid_radius_mean, row$selected_grid_radius_mcse, 3L
      ), " & ",
      panic_cell(row$mean_shared_runtime, row$mcse_shared_runtime, 4L),
      " \\\\"
    ))
  }
  c(lines, "\\hline", "\\end{tabular}", "\\end{adjustbox}",
    "\\end{table}")
}

render_boundary_table <- function(boundary_summary) {
  required <- c(
    "n", "replications", "empirical_atom", "atom_mcse",
    "theoretical_atom", "q10", "q25", "theoretical_q10",
    "theoretical_q25"
  )
  panic_require_columns(boundary_summary, required, "boundary_summary")
  rows <- boundary_summary[order(boundary_summary$n), , drop = FALSE]
  replications <- panic_single_value(
    rows$replications, "boundary replications"
  )
  theoretical_atom <- panic_single_value(
    rows$theoretical_atom, "boundary theoretical atom"
  )
  theoretical_q10 <- panic_single_value(
    rows$theoretical_q10, "boundary theoretical tenth percentile"
  )
  theoretical_q25 <- panic_single_value(
    rows$theoretical_q25, "boundary theoretical twenty-fifth percentile"
  )

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0(
      "\\caption{Scalar binding-boundary illustration with ",
      panic_plain_integer(replications), " replications. The limiting law has ",
      "a point mass $1/2$ at zero and a standard-Gaussian density on the ",
      "negative half-line. MCSE denotes the binomial Monte Carlo standard error ",
      "of the atom estimate.}"
    ),
    "\\label{Table: boundary illustration}",
    "\\small",
    "\\begin{tabular}{rrrrrr}",
    "\\hline",
    "$n$ & Atom & MCSE & 10th percentile & 25th percentile & Limit atom \\\\",
    "\\hline"
  )
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, , drop = FALSE]
    lines <- c(lines, paste0(
      sprintf(
        "%d & %.4f & %.4f & %.3f & %.3f & %.1f",
        as.integer(row$n), row$empirical_atom, row$atom_mcse, row$q10,
        row$q25, row$theoretical_atom
      ),
      " \\\\"
    ))
  }
  limit_line <- paste0(
    sprintf(
      "Limit & %.1f & -- & %.3f & %.3f & %.1f",
      theoretical_atom, theoretical_q10, theoretical_q25, theoretical_atom
    ),
    " \\\\"
  )
  c(
    lines,
    limit_line,
    "\\hline", "\\end{tabular}", "\\end{table}"
  )
}

render_subsampling_table <- function(subsampling_summary, boundary_summary) {
  required <- c(
    "n", "b", "b_over_n", "random_subsamples", "atom_at_zero", "q10",
    "q25", "cdf_minus_1", "cdf_minus_0_5", "cdf_minus_0_1",
    "theoretical_cdf_minus_1", "theoretical_cdf_minus_0_5",
    "theoretical_cdf_minus_0_1"
  )
  panic_require_columns(subsampling_summary, required, "subsampling_summary")
  panic_require_columns(
    boundary_summary,
    c("theoretical_atom", "theoretical_q10", "theoretical_q25"),
    "boundary_summary"
  )
  rows <- subsampling_summary[order(subsampling_summary$b), , drop = FALSE]
  n <- panic_single_value(rows$n, "subsampling full-sample size")
  replications <- panic_single_value(
    rows$random_subsamples, "random subsamples"
  )
  theoretical_atom <- panic_single_value(
    boundary_summary$theoretical_atom, "boundary theoretical atom"
  )
  theoretical_q10 <- panic_single_value(
    boundary_summary$theoretical_q10, "boundary theoretical tenth percentile"
  )
  theoretical_q25 <- panic_single_value(
    boundary_summary$theoretical_q25,
    "boundary theoretical twenty-fifth percentile"
  )
  theoretical_cdf_minus_1 <- panic_single_value(
    rows$theoretical_cdf_minus_1, "theoretical cdf at -1"
  )
  theoretical_cdf_minus_0_5 <- panic_single_value(
    rows$theoretical_cdf_minus_0_5, "theoretical cdf at -0.5"
  )
  theoretical_cdf_minus_0_1 <- panic_single_value(
    rows$theoretical_cdf_minus_0_1, "theoretical cdf at -0.1"
  )

  lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0(
      "\\caption{Subsampling sensitivity for one binding full sample of size ",
      "$n=", sprintf("%d", as.integer(n)), "$, using ",
      panic_plain_integer(replications), " random subsamples at each $b$. The ",
      "final three columns compare empirical cdf values with the limit cdf at ",
      "continuity points.}"
    ),
    "\\label{Table: subsampling sensitivity}",
    "\\begin{adjustbox}{width=\\textwidth}",
    "\\small",
    "\\begin{tabular}{rrrrrrrr}",
    "\\hline",
    paste0(
      "$b$ & $b/n$ & Atom at zero & 10th pct. & 25th pct. & $F(-1)$ & ",
      "$F(-0.5)$ & $F(-0.1)$ \\\\"
    ),
    "\\hline"
  )
  for (i in seq_len(nrow(rows))) {
    row <- rows[i, , drop = FALSE]
    lines <- c(lines, paste0(
      sprintf(
        "%d & %.2f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f",
        as.integer(row$b), row$b_over_n, row$atom_at_zero, row$q10, row$q25,
        row$cdf_minus_1, row$cdf_minus_0_5, row$cdf_minus_0_1
      ),
      " \\\\"
    ))
  }
  limit_line <- paste0(
    sprintf(
      paste0(
        "Limit & 0 & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f"
      ),
      theoretical_atom, theoretical_q10, theoretical_q25,
      theoretical_cdf_minus_1, theoretical_cdf_minus_0_5,
      theoretical_cdf_minus_0_1
    ),
    " \\\\"
  )
  c(
    lines,
    limit_line,
    "\\hline", "\\end{tabular}", "\\end{adjustbox}", "\\end{table}"
  )
}

panic_read_result_csv <- function(results_dir, filename) {
  path <- file.path(results_dir, filename)
  if (!file.exists(path)) {
    stop("Missing required result summary: ", path, call. = FALSE)
  }
  read.csv(
    path, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN")
  )
}

render_all_manuscript_tables <- function(results_dir) {
  configuration <- panic_read_result_csv(results_dir, "configuration.csv")
  simulation <- panic_read_result_csv(results_dir, "simulation_summary.csv")
  decision <- panic_read_result_csv(results_dir, "confirmatory_decision.csv")
  targets <- panic_read_result_csv(
    results_dir, "calibration_target_summary.csv"
  )
  diagnostics <- panic_read_result_csv(results_dir, "diagnostic_summary.csv")
  grid <- panic_read_result_csv(results_dir, "grid_sensitivity_summary.csv")
  boundary <- panic_read_result_csv(results_dir, "boundary_summary.csv")
  subsampling <- panic_read_result_csv(results_dir, "subsampling_summary.csv")

  setNames(list(
    render_primary_support_table(simulation, configuration),
    render_primary_performance_table(simulation, configuration),
    render_confirmatory_decision_table(decision),
    render_calibration_table(targets, diagnostics),
    render_grid_sensitivity_table(grid),
    render_boundary_table(boundary),
    render_subsampling_table(subsampling, boundary)
  ), PANIC_RENDERED_TABLE_FILES)
}

panic_table_bytes <- function(lines) {
  if (!is.character(lines) || !length(lines) || anyNA(lines) ||
      any(grepl("\r|\n", lines))) {
    stop("Rendered table must be a non-empty vector of newline-free strings",
         call. = FALSE)
  }
  charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n")))
}

write_manuscript_table_subset <- function(tables, output_dir) {
  if (!is.list(tables) || !length(tables) || is.null(names(tables)) ||
      any(!nzchar(names(tables))) || anyDuplicated(names(tables)) ||
      !all(names(tables) %in% PANIC_RENDERED_TABLE_FILES)) {
    stop(
      "Rendered table subset must be a non-empty, uniquely named subset of ",
      "PANIC_RENDERED_TABLE_FILES", call. = FALSE
    )
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }
  output_paths <- file.path(output_dir, names(tables))
  for (i in seq_along(tables)) {
    connection <- file(output_paths[[i]], open = "wb")
    on.exit(close(connection), add = TRUE)
    writeBin(panic_table_bytes(tables[[i]]), connection)
    close(connection)
    on.exit(NULL, add = FALSE)
  }
  invisible(output_paths)
}

write_all_manuscript_tables <- function(tables, output_dir) {
  if (!identical(names(tables), PANIC_RENDERED_TABLE_FILES)) {
    stop("Rendered table list has unexpected names or order", call. = FALSE)
  }
  write_manuscript_table_subset(tables, output_dir)
}
