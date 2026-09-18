#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]))
} else normalizePath("Simulation_Study.R")
production_dir <- dirname(script_path)

run <- function(script, extra = character()) {
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(file.path(production_dir, script), extra)
  )
  if (status != 0L) stop(script, " failed with status ", status)
}

cat("1. Validate the PanIC-CF implementation.\n")
run("Implementation_Validation.R", args[grepl("^--output-dir=", args)])
cat("2. Run the confirmatory simulations.\n")
run("Run_Confirmatory_Simulations.R", args)
cat("3. Generate manuscript summaries.\n")
run("Simulation_Results_Analysis.R", args[grepl("^--output-dir=", args)])
cat("4. Verify the generated replication-level results.\n")
run("Results_Verification.R", args[grepl("^--output-dir=", args)])
