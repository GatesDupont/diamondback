# Durable results on disk: a directory with the labelled raster, the metric
# table, and the metadata that identifies the run. This is also the cache
# format, so a cached result and a saved result are the same thing.

RESULT_FILES <- list(
  labels = "patches.tif",
  metrics = "metrics.csv",
  metadata = "metadata.json"
)

#' Write a patch result to disk
#'
#' Saves a [patch_result] as a small, durable directory: the labelled raster as
#' GeoTIFF, the metric table as CSV, and the metadata as JSON. Nothing is
#' pickled and nothing depends on an R or Python session, so results stay
#' readable years later and from other tools.
#'
#' This is also the cache format used by [analyze_patches()], so a result you
#' saved and a result it cached are interchangeable.
#'
#' @param x A [patch_result].
#' @param dir Directory to write to; created if needed.
#' @param overwrite Overwrite existing files in `dir`.
#' @param quiet Suppress progress reporting.
#' @return `dir`, invisibly.
#' @seealso [read_patch_result()]
#' @export
#' @examplesIf diamondback_ready()
#' m <- matrix(c(1, 1, 0, 0, 1, 0, 0, 1, 1), nrow = 3)
#' res <- analyze_patches(m, quiet = TRUE)
#' d <- file.path(tempdir(), "patches_demo")
#' write_patch_result(res, d, overwrite = TRUE, quiet = TRUE)
#' read_patch_result(d)
write_patch_result <- function(x, dir, overwrite = FALSE, quiet = FALSE) {
  stopifnot(inherits(x, "patch_result"))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  f_lab <- file.path(dir, RESULT_FILES$labels)
  f_met <- file.path(dir, RESULT_FILES$metrics)
  f_meta <- file.path(dir, RESULT_FILES$metadata)

  existing <- c(f_lab, f_met, f_meta)
  existing <- existing[file.exists(existing)]
  if (length(existing) && !isTRUE(overwrite)) {
    cli::cli_abort(c(
      "{.path {dir}} already contains a patch result.",
      "x" = "Found {.file {basename(existing)}}.",
      "i" = "Use {.code overwrite = TRUE} to replace it."
    ), call = NULL)
  }

  p <- .subset2(x, "path")
  if (!is.null(p) && normalizePath(p, mustWork = FALSE) == normalizePath(f_lab, mustWork = FALSE)) {
    # Already written here by label_patches(output = ...); nothing to copy.
  } else {
    r <- db_patches(x)
    db_register_file(f_lab)
    terra::writeRaster(r, f_lab, overwrite = TRUE, datatype = "INT4S",
                       NAflag = -2147483648)
  }

  utils::write.csv(.subset2(x, "metrics"), f_met, row.names = FALSE)

  meta <- .subset2(x, "metadata")
  meta$files <- RESULT_FILES
  # Store the canonical cache fields, not just the key: they survive JSON
  # exactly, and keeping them means a cache miss can still name the field that
  # differed instead of just saying "something changed".
  meta$cache_fields <- db_cache_fields(meta)
  meta$cache_key <- db_cache_key(meta)
  # digits = NA writes full double precision. Rounding here would change a
  # lon/lat resolution enough to break the cache key on the way back in.
  writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, null = "null",
                              na = "string", digits = NA, pretty = TRUE), f_meta)

  if (!quiet) cli::cli_alert_success("Wrote patch result to {.path {dir}}.")
  invisible(dir)
}

#' Read a patch result from disk
#'
#' Reads back what [write_patch_result()] wrote. The labelled raster is opened
#' from the file rather than loaded, so reading a result is cheap regardless of
#' its size.
#'
#' The NumPy arrays from the original run are gone, so the first call to
#' [patch_metrics()] or [compare_patches()] on a restored result rebuilds them
#' from the raster. Everything works; the first such call just pays a read.
#'
#' @param dir Directory written by [write_patch_result()].
#' @return A [patch_result].
#' @seealso [write_patch_result()]
#' @export
#' @examplesIf diamondback_ready()
#' m <- matrix(c(1, 1, 0, 0, 1, 0, 0, 1, 1), nrow = 3)
#' d <- file.path(tempdir(), "patches_demo2")
#' write_patch_result(analyze_patches(m, quiet = TRUE), d, overwrite = TRUE, quiet = TRUE)
#' read_patch_result(d)
read_patch_result <- function(dir) {
  f_lab <- file.path(dir, RESULT_FILES$labels)
  f_met <- file.path(dir, RESULT_FILES$metrics)
  f_meta <- file.path(dir, RESULT_FILES$metadata)

  missing <- c(f_lab, f_met, f_meta)[!file.exists(c(f_lab, f_met, f_meta))]
  if (length(missing)) {
    cli::cli_abort(c(
      "{.path {dir}} is not a complete patch result.",
      "x" = "Missing {.file {basename(missing)}}.",
      "i" = "Results are written by {.fn write_patch_result} or by \\
             {.code analyze_patches(output_dir = ...)}."
    ), call = NULL)
  }

  meta <- jsonlite::fromJSON(f_meta, simplifyVector = TRUE)
  metrics <- utils::read.csv(f_met)

  new_patch_result(
    patches = NULL,
    patches_path = f_lab,
    metrics = metrics,
    metadata = meta,
    arrays = NULL
  )
}
