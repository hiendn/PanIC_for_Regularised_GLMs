#!/usr/bin/env Rscript

## Refresh and verify the self-contained manuscript table blocks.
##
## The files under results/table_*.tex remain the auditable, generated table
## fragments.  The manuscript embeds exact copies between stable marker lines,
## so LaTeX compilation has no external table-file dependency while the code
## can still detect any drift between a result and its displayed table.

PANIC_TABLE_FILES <- c(
  "table_primary_support.tex",
  "table_primary_performance.tex",
  "table_calibration.tex",
  "table_grid_sensitivity.tex",
  "table_boundary.tex",
  "table_subsampling.tex"
)

panic_table_marker <- function(kind, name) {
  sprintf("%% %s PANIC INLINE TABLE: %s", kind, name)
}

panic_marker_inventory <- function(manuscript_lines) {
  pattern <- "^% (BEGIN|END) PANIC INLINE TABLE: (.+)$"
  marker_lines <- grep(pattern, manuscript_lines, value = TRUE)
  matches <- regmatches(marker_lines, regexec(pattern, marker_lines))
  kinds <- if (length(matches)) {
    vapply(matches, `[[`, character(1), 2L)
  } else {
    character()
  }
  names <- if (length(matches)) {
    vapply(matches, `[[`, character(1), 3L)
  } else {
    character()
  }
  begin_names <- names[kinds == "BEGIN"]
  end_names <- names[kinds == "END"]
  expected <- sort(PANIC_TABLE_FILES)
  passed <-
    length(begin_names) == length(PANIC_TABLE_FILES) &&
    length(end_names) == length(PANIC_TABLE_FILES) &&
    !anyDuplicated(begin_names) && !anyDuplicated(end_names) &&
    identical(sort(begin_names), expected) &&
    identical(sort(end_names), expected)
  list(
    passed = passed,
    detail = paste0(
      "BEGIN={", paste(begin_names, collapse = ","), "};END={",
      paste(end_names, collapse = ","), "}"
    )
  )
}

read_text_lines <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

read_binary_bytes <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  size <- file.info(path)$size
  readBin(path, what = "raw", n = size)
}

locate_raw_sequence <- function(bytes, sequence, description) {
  if (!length(sequence) || length(sequence) > length(bytes)) {
    stop(description, " must occur exactly once; found 0", call. = FALSE)
  }
  candidates <- which(bytes == sequence[[1L]])
  hits <- candidates[vapply(candidates, function(first) {
    last <- first + length(sequence) - 1L
    last <= length(bytes) && identical(bytes[first:last], sequence)
  }, logical(1))]
  if (length(hits) != 1L) {
    stop(
      description, " must occur exactly once; found ", length(hits),
      call. = FALSE
    )
  }
  hits[[1L]]
}

locate_unique_line <- function(lines, value, description) {
  index <- which(lines == value)
  if (length(index) != 1L) {
    stop(
      description, " must occur exactly once; found ", length(index),
      call. = FALSE
    )
  }
  index
}

validate_table_fragment <- function(lines, name) {
  if (!length(lines) || lines[[1L]] != "\\begin{table}[H]" ||
      lines[[length(lines)]] != "\\end{table}") {
    stop(name, " is not a complete [H] table float", call. = FALSE)
  }
  if (any(grepl("TODO", lines, ignore.case = TRUE))) {
    stop(name, " contains a TODO marker", call. = FALSE)
  }
  invisible(TRUE)
}

extract_inline_table <- function(manuscript_lines, name) {
  begin <- locate_unique_line(
    manuscript_lines, panic_table_marker("BEGIN", name),
    paste0("BEGIN marker for ", name)
  )
  end <- locate_unique_line(
    manuscript_lines, panic_table_marker("END", name),
    paste0("END marker for ", name)
  )
  if (end <= begin + 1L) {
    stop("Empty or reversed inline block for ", name, call. = FALSE)
  }
  manuscript_lines[seq.int(begin + 1L, end - 1L)]
}

extract_inline_table_bytes <- function(manuscript_path, name) {
  manuscript_bytes <- read_binary_bytes(manuscript_path)
  begin_sequence <- charToRaw(paste0(
    panic_table_marker("BEGIN", name), "\n"
  ))
  end_sequence <- charToRaw(panic_table_marker("END", name))
  begin <- locate_raw_sequence(
    manuscript_bytes, begin_sequence, paste0("BEGIN marker for ", name)
  )
  end <- locate_raw_sequence(
    manuscript_bytes, end_sequence, paste0("END marker for ", name)
  )
  first <- begin + length(begin_sequence)
  last <- end - 1L
  if (last < first) {
    stop("Empty or reversed inline byte block for ", name, call. = FALSE)
  }
  manuscript_bytes[first:last]
}

replace_line_range <- function(lines, first, last, replacement) {
  before <- if (first > 1L) lines[seq_len(first - 1L)] else character()
  after <- if (last < length(lines)) lines[seq.int(last + 1L, length(lines))]
    else character()
  c(before, replacement, after)
}

refresh_one_inline_table <- function(manuscript_lines, table_lines, name) {
  validate_table_fragment(table_lines, name)
  begin_marker <- panic_table_marker("BEGIN", name)
  end_marker <- panic_table_marker("END", name)
  replacement <- c(begin_marker, table_lines, end_marker)

  begin_hits <- which(manuscript_lines == begin_marker)
  end_hits <- which(manuscript_lines == end_marker)
  if (length(begin_hits) || length(end_hits)) {
    if (length(begin_hits) != 1L || length(end_hits) != 1L ||
        end_hits <= begin_hits) {
      stop("Malformed existing inline markers for ", name, call. = FALSE)
    }
    return(replace_line_range(
      manuscript_lines, begin_hits, end_hits, replacement
    ))
  }

  ## Migration path for the audited working manuscript, whose table inputs are
  ## wrapped in IfFileExists fallbacks.
  legacy_start <- sprintf(
    "\\IfFileExists{%s}{\\input{%s}}{%%", name, name
  )
  start_hits <- which(manuscript_lines == legacy_start)
  if (length(start_hits) == 1L) {
    candidates <- which(
      seq_along(manuscript_lines) >= start_hits &
        manuscript_lines == "\\end{table}}"
    )
    if (!length(candidates)) {
      stop("Cannot find the end of the legacy fallback for ", name,
           call. = FALSE)
    }
    return(replace_line_range(
      manuscript_lines, start_hits, candidates[[1L]], replacement
    ))
  }
  if (length(start_hits) > 1L) {
    stop("Legacy fallback occurs more than once for ", name, call. = FALSE)
  }

  ## Migration path for the previous manuscript, which used a plain input.
  plain_input <- sprintf("\\input{%s}", name)
  input_hits <- which(trimws(manuscript_lines) == plain_input)
  if (length(input_hits) == 1L) {
    return(replace_line_range(
      manuscript_lines, input_hits, input_hits, replacement
    ))
  }

  stop("No replaceable table location found for ", name, call. = FALSE)
}

verify_inline_tables <- function(manuscript_path, tables_dir,
                                 stop_on_failure = TRUE) {
  manuscript_lines <- read_text_lines(manuscript_path)
  marker_inventory <- panic_marker_inventory(manuscript_lines)
  details <- character(length(PANIC_TABLE_FILES))
  passed <- logical(length(PANIC_TABLE_FILES))

  for (i in seq_along(PANIC_TABLE_FILES)) {
    name <- PANIC_TABLE_FILES[[i]]
    expected <- read_text_lines(file.path(tables_dir, name))
    expected_bytes <- read_binary_bytes(file.path(tables_dir, name))
    validate_table_fragment(expected, name)
    observed <- tryCatch(
      extract_inline_table_bytes(manuscript_path, name),
      error = function(error) error
    )
    passed[[i]] <- !inherits(observed, "error") &&
      identical(observed, expected_bytes)
    details[[i]] <- if (inherits(observed, "error")) {
      conditionMessage(observed)
    } else if (!identical(observed, expected_bytes)) {
      "embedded bytes differ from the retained result fragment"
    } else {
      paste0(length(expected), " lines and ", length(expected_bytes),
             " bytes match exactly")
    }
  }

  table_reference_pattern <- paste0(
    "\\\\(?:input|include)\\s*\\{[^}]*table_[^}]*\\}",
    "|\\\\IfFileExists\\s*\\{table_[^}]*\\}"
  )
  no_external_references <- !any(grepl(
    table_reference_pattern, manuscript_lines, perl = TRUE
  ))
  result <- data.frame(
    table = PANIC_TABLE_FILES,
    passed = as.integer(passed),
    detail = details,
    stringsAsFactors = FALSE
  )
  attr(result, "no_external_table_references") <- no_external_references
  attr(result, "exact_marker_inventory") <- marker_inventory$passed
  attr(result, "marker_inventory_detail") <- marker_inventory$detail

  if (isTRUE(stop_on_failure) &&
      (!all(passed) || !no_external_references ||
       !marker_inventory$passed)) {
    failed <- result$table[result$passed != 1L]
    extra <- if (!no_external_references) "external table reference remains"
      else character()
    marker_problem <- if (!marker_inventory$passed) {
      paste0("marker inventory is not exactly the six-table contract: ",
             marker_inventory$detail)
    } else {
      character()
    }
    stop(
      "Inline-table verification failed: ",
      paste(c(failed, extra, marker_problem), collapse = "; "),
      call. = FALSE
    )
  }
  result
}

refresh_inline_tables <- function(manuscript_path, tables_dir) {
  manuscript_lines <- read_text_lines(manuscript_path)
  for (name in PANIC_TABLE_FILES) {
    table_lines <- read_text_lines(file.path(tables_dir, name))
    manuscript_lines <- refresh_one_inline_table(
      manuscript_lines, table_lines, name
    )
  }

  temporary <- tempfile(
    pattern = paste0(basename(manuscript_path), "."),
    tmpdir = dirname(manuscript_path)
  )
  on.exit(unlink(temporary), add = TRUE)
  writeLines(manuscript_lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, manuscript_path)) {
    stop("Could not atomically replace ", manuscript_path, call. = FALSE)
  }
  verify_inline_tables(manuscript_path, tables_dir, stop_on_failure = TRUE)
}

parse_named_argument <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", prefix), "", hit[[1L]])
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if ("--refresh" %in% args) "refresh" else "verify"
  manuscript_path <- parse_named_argument(args, "manuscript")
  tables_dir <- parse_named_argument(args, "tables-dir")
  if (is.null(manuscript_path) || is.null(tables_dir)) {
    stop(
      "Usage: Rscript Manuscript_Table_Tools.R [--refresh|--verify] ",
      "--manuscript=PATH --tables-dir=PATH", call. = FALSE
    )
  }
  if (mode == "refresh") {
    result <- refresh_inline_tables(manuscript_path, tables_dir)
    cat("Refreshed and verified six inline manuscript tables.\n")
  } else {
    result <- verify_inline_tables(
      manuscript_path, tables_dir, stop_on_failure = TRUE
    )
    cat("Verified six inline manuscript tables.\n")
  }
  print(result, row.names = FALSE)
}
