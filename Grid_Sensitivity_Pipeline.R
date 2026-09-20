#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]))
} else {
  normalizePath("Grid_Sensitivity_Pipeline.R")
}
candidate_dir <- dirname(script_path)

run <- function(script, extra = character()) {
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      shQuote(file.path(candidate_dir, script)),
      vapply(extra, shQuote, character(1))
    )
  )
  if (status != 0L) stop(script, " failed with status ", status)
}

output_arg <- args[grepl("^--output-dir=", args)]
run("Implementation_Validation.R")
run("Grid_Sensitivity_Study.R", args)
run("Grid_Results_Verification.R", output_arg)
