## Verify the immutable source manifest before a simulation is executed.

verify_method_lock <- function(candidate_dir) {
  manifest_path <- file.path(candidate_dir, "METHOD_LOCK_SHA256.txt")
  if (!file.exists(manifest_path)) {
    stop("Missing METHOD_LOCK_SHA256.txt", call. = FALSE)
  }
  lines <- readLines(manifest_path, warn = FALSE)
  entries <- grep("^[0-9a-f]{64}  ", lines, value = TRUE)
  if (!length(entries)) stop("The method-lock manifest is empty", call. = FALSE)
  expected <- substr(entries, 1L, 64L)
  relative <- substring(entries, 67L)
  paths <- normalizePath(
    file.path(candidate_dir, relative), mustWork = TRUE
  )
  actual <- vapply(paths, function(path) {
    output <- system2(
      "shasum", c("-a", "256", shQuote(path)),
      stdout = TRUE, stderr = TRUE
    )
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) {
      stop("shasum failed for ", path, call. = FALSE)
    }
    substr(output[[1L]], 1L, 64L)
  }, character(1))
  mismatch <- which(actual != expected)
  if (length(mismatch)) {
    stop(
      "Method-lock mismatch: ", paste(relative[mismatch], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(data.frame(file = relative, sha256 = actual))
}
