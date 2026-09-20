#!/usr/bin/env Rscript

script_argument <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
script_path <- if (length(script_argument)) {
  sub("^--file=", "", script_argument[[1L]])
} else {
  "Render_Manuscript_Tables.R"
}
script_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
source(file.path(script_dir, "Manuscript_Table_Rendering.R"))

parse_cli_arguments <- function(arguments) {
  values <- list()
  index <- 1L
  while (index <= length(arguments)) {
    argument <- arguments[[index]]
    if (argument %in% c("--help", "-h")) {
      values$help <- TRUE
      index <- index + 1L
      next
    }
    matched <- regexec("^--(results-dir|output-dir)=(.*)$", argument)
    pieces <- regmatches(argument, matched)[[1L]]
    if (length(pieces)) {
      values[[pieces[[2L]]]] <- pieces[[3L]]
      index <- index + 1L
      next
    }
    if (argument %in% c("--results-dir", "--output-dir")) {
      if (index == length(arguments)) {
        stop("Missing value after ", argument, call. = FALSE)
      }
      values[[sub("^--", "", argument)]] <- arguments[[index + 1L]]
      index <- index + 2L
      next
    }
    stop("Unknown argument: ", argument, call. = FALSE)
  }
  values
}

usage <- paste0(
  "Usage: Rscript Render_Manuscript_Tables.R ",
  "--results-dir PATH --output-dir PATH\n"
)

arguments <- parse_cli_arguments(commandArgs(trailingOnly = TRUE))
if (isTRUE(arguments$help)) {
  cat(usage)
  quit(save = "no", status = 0L)
}
if (is.null(arguments[["results-dir"]]) ||
    is.null(arguments[["output-dir"]]) ||
    !nzchar(arguments[["results-dir"]]) ||
    !nzchar(arguments[["output-dir"]])) {
  stop(usage, call. = FALSE)
}

results_dir <- normalizePath(
  arguments[["results-dir"]], mustWork = TRUE
)
output_dir <- arguments[["output-dir"]]
tables <- render_all_manuscript_tables(results_dir)
paths <- write_all_manuscript_tables(tables, output_dir)

cat(
  "Rendered ", length(paths), " manuscript table fragments to ",
  normalizePath(output_dir, mustWork = TRUE), ".\n", sep = ""
)
