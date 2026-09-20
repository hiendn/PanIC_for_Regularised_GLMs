## Render the revised runtime figure from archived replication-level timings.

if (!exists("production_dir", inherits = FALSE)) {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg)) {
    normalizePath(sub("^--file=", "", script_arg[[1L]]))
  } else normalizePath("Render_Runtime_Figure.R")
  production_dir <- dirname(script_path)
}
if (!exists("results_dir", inherits = FALSE)) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep("^--output-dir=", args, value = TRUE)
  output_name <- if (length(hit)) sub("^--output-dir=", "", hit[[1L]]) else
    "results"
  results_dir <- if (grepl("^/", output_name)) output_name else
    file.path(production_dir, output_name)
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
runtime <- read.csv(file.path(results_dir, "runtime_raw.csv"),
                    stringsAsFactors = FALSE)
families <- c("gaussian", "binomial", "poisson")

pdf(
  file.path(results_dir, "figure_runtime.pdf"),
  width = 7.2, height = 3.8, family = "Helvetica", pointsize = 9
)
par(mfrow = c(1, 3), mar = c(5.5, 3.8, 2.0, 0.7),
    mgp = c(2.3, 0.7, 0), tcl = -0.25)
labels <- c(gaussian = "Gaussian", binomial = "Logistic", poisson = "Poisson")
cols <- c("#2166AC", "#67A9CF", "#B2182B")
for (family in families) {
  subset <- runtime[runtime$family == family & runtime$failed == 0L, ]
  methods <- c(
    "Revised PanIC-CF", "5-fold CV",
    if (family == "gaussian") "BIC-like" else "Exploratory active-count"
  )
  split_values <- lapply(methods, function(method) {
    subset$elapsed_seconds[subset$method == method]
  })
  method_labels <- c(
    "PanIC-CF", "Shared\n5-fold CV", "BIC-like"
  )
  boxplot(
    split_values, names = method_labels,
    col = cols, border = "grey30", outline = TRUE,
    outpch = 1, outcex = 0.55, las = 2, cex.axis = 0.82,
    ylab = if (family == families[1L]) "elapsed seconds" else "",
    main = labels[[family]]
  )
}
dev.off()
