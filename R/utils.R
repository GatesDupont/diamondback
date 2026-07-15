# Progress reporting, timing, fingerprints, and the temp-file registry.

# ---------------------------------------------------------------------------
# progress
# ---------------------------------------------------------------------------

# Stage-level progress only. Per-cell progress would flood the console and cost
# more than it tells you; what a user waiting on a long run needs to know is
# which stage they are in and how long it has taken.

db_progress <- function(name, total, quiet = FALSE) {
  if (isTRUE(quiet) || total <= 1L) return(NULL)
  cli::cli_progress_bar(name, total = total, .envir = parent.frame(2))
}

db_progress_step <- function(pb, i) {
  if (!is.null(pb)) cli::cli_progress_update(id = pb, set = i)
  invisible(NULL)
}

db_progress_done <- function(pb) {
  if (!is.null(pb)) try(cli::cli_progress_done(id = pb), silent = TRUE)
  invisible(NULL)
}

db_elapsed_num <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

db_elapsed <- function(t0) {
  s <- db_elapsed_num(t0)
  if (s < 1) return(sprintf("%.0f ms", s * 1000))
  if (s < 60) return(sprintf("%.1f s", s))
  if (s < 3600) return(sprintf("%.1f min", s / 60))
  sprintf("%.2f h", s / 3600)
}

# ---------------------------------------------------------------------------
# fingerprints
# ---------------------------------------------------------------------------

# How large a file may be before "auto" stops hashing its contents and falls
# back to size + mtime. Hashing a 40 GB raster to decide whether to reuse a
# cache can cost more than the labelling did.
HASH_MAX_MB <- 200

#' Describe an input well enough to detect that it has changed
#'
#' `fingerprint` controls the strength/cost trade-off:
#'
#' * `"auto"` (default) --- content hash for files under `hash_max_mb`, size and
#'   mtime beyond that.
#' * `"full"` --- always hash the contents, however large. Slower, but immune to
#'   a file being rewritten with the same size and a restored mtime.
#' * `"fast"` --- never hash; size and mtime only.
#'
#' Size and mtime are a *heuristic*: a file can be rewritten at the same size
#' with `touch -r` restoring its timestamp, and the cache would then be reused
#' for changed data. That is why `"full"` exists, and why the mode is recorded
#' in the metadata and forms part of the cache key -- a cache built under a weak
#' fingerprint is never silently trusted by a request asking for a strong one.
#'
#' For in-memory rasters there is no file to fingerprint. Small ones are hashed
#' by value; large ones are marked unfingerprintable, and caching is declined
#' rather than guessed at.
#' @noRd
db_source_info <- function(x, fingerprint = c("auto", "full", "fast"),
                           hash_max_mb = HASH_MAX_MB) {
  fingerprint <- match.arg(fingerprint)

  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    fi <- file.info(x)
    lim <- switch(fingerprint, auto = hash_max_mb, full = Inf, fast = -1)
    return(list(
      type = "file",
      path = normalizePath(x, winslash = "/"),
      size = as.numeric(fi$size),
      mtime = as.character(fi$mtime),
      hash = db_file_hash(x, lim),
      fingerprint = fingerprint
    ))
  }

  if (inherits(x, "SpatRaster")) {
    srcs <- terra::sources(x)
    srcs <- srcs[nzchar(srcs)]
    if (length(srcs) == 1L && file.exists(srcs)) {
      return(db_source_info(srcs, fingerprint, hash_max_mb))
    }
    # An in-memory raster has no path, size or mtime to fall back on, so the
    # only identity available is a hash of its values. That makes the modes
    # mean something slightly different here than for a file:
    #   full  - always hash, however large
    #   auto  - hash while it is cheap (<= 1e6 cells), otherwise no identity
    #   fast  - never hash, so no identity at all
    # "No identity" means caching is declined rather than guessed at; see
    # db_cache_match(). This is the safe direction: terra drops a raster's file
    # source as soon as it is modified in memory, so anything reaching here is
    # data we have never seen written down.
    lim <- switch(fingerprint, full = Inf, auto = 1e6, fast = -1)
    if (terra::ncell(x) <= lim) {
      return(list(
        type = "memory",
        path = NA_character_,
        size = terra::ncell(x),
        mtime = NA_character_,
        hash = digest::digest(terra::values(x, mat = FALSE), algo = "sha1"),
        fingerprint = fingerprint
      ))
    }
    return(list(type = "memory", path = NA_character_, size = terra::ncell(x),
                mtime = NA_character_, hash = NA_character_,
                fingerprint = fingerprint))
  }

  if (is.matrix(x)) {
    return(list(type = "matrix", path = NA_character_, size = length(x),
                mtime = NA_character_, hash = digest::digest(x, algo = "sha1"),
                fingerprint = fingerprint))
  }

  list(type = "unknown", path = NA_character_, size = NA_real_,
       mtime = NA_character_, hash = NA_character_, fingerprint = fingerprint)
}

db_file_hash <- function(path, hash_max_mb = HASH_MAX_MB) {
  if (is.infinite(hash_max_mb) && hash_max_mb < 0) return(NA_character_)
  size_mb <- file.info(path)$size / 1024^2
  if (is.na(size_mb) || size_mb > hash_max_mb) return(NA_character_)
  tryCatch(digest::digest(file = path, algo = "sha1"), error = function(e) NA_character_)
}

#' Assemble the metadata that both documents and identifies a run
#' @noRd
db_metadata <- function(source, mask_source, geometry, class, directions, na,
                        crop, counts, n_patches, elapsed) {
  vers <- tryCatch(py_r(db_py()$versions()), error = function(e) list())
  list(
    package_version = as.character(utils::packageVersion("diamondback")),
    algorithm_version = ALGORITHM_VERSION,
    created = as.character(Sys.time()),
    elapsed_secs = elapsed,
    source = source,
    mask_source = mask_source,
    geometry = geometry,
    class = if (is.null(class)) NA_real_ else class,
    binary = is.null(class),
    directions = directions,
    na = na,
    crop = crop,
    n_patches = n_patches,
    cells = list(
      total = geometry$ncell,
      foreground = counts$foreground,
      background = counts$background,
      missing = counts$missing,
      outside = counts$outside
    ),
    backend = list(
      python = vers$python %||% NA_character_,
      numpy = vers$numpy %||% NA_character_,
      scipy = vers$scipy %||% NA_character_
    ),
    validated = NA
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# temp-file registry
# ---------------------------------------------------------------------------

# Files diamondback created in this session. Cleanup only ever touches this
# list. terra::tmpFiles(remove = TRUE) is never called: it deletes files backing
# SpatRaster objects that other code may still be holding, which is exactly the
# failure the old pipelines kept hitting.
.db_files <- new.env(parent = emptyenv())
.db_files$created <- character()   # named: value = path, name = "temp"/"output"

#' Record a file diamondback created
#'
#' `kind` is recorded, not inferred from the path. A user's `output` can
#' perfectly well live under `tempdir()` -- that makes it a temporary location,
#' not a file we are entitled to delete.
#' @noRd
db_register_file <- function(path, kind = c("output", "temp")) {
  kind <- match.arg(kind)
  keep <- .db_files$created[.db_files$created != path]
  .db_files$created <- c(keep, stats::setNames(path, kind))
  invisible(path)
}

#' Make a temporary raster file that diamondback owns
#' @noRd
db_tempfile <- function(ext = ".tif") {
  f <- tempfile(pattern = "diamondback_", fileext = ext)
  db_register_file(f, kind = "temp")
  f
}

#' Remove temporary files created by diamondback in this session
#'
#' Only files this package created are removed, and only those still on disk.
#' Files you passed as `output` are yours: they are tracked so that they can be
#' reported, but they are not deleted unless `include_outputs = TRUE`.
#'
#' This exists because the usual reflex --- `terra::tmpFiles(remove = TRUE)` ---
#' deletes files backing `SpatRaster` objects that are still live, which
#' corrupts unrelated work. diamondback never calls it.
#'
#' @param include_outputs Also remove files written via an `output` argument.
#' @param quiet Suppress the report.
#' @return Character vector of removed paths, invisibly.
#' @export
#' @examples
#' db_clean_temp()
db_clean_temp <- function(include_outputs = FALSE, quiet = FALSE) {
  files <- .db_files$created
  if (!isTRUE(include_outputs)) {
    files <- files[names(files) == "temp"]
  }
  files <- files[file.exists(files)]
  if (!length(files)) {
    if (!quiet) cli::cli_alert_info("No diamondback temporary files to remove.")
    return(invisible(character()))
  }
  ok <- file.remove(files)
  removed <- unname(files[ok])
  .db_files$created <- .db_files$created[!(.db_files$created %in% removed)]
  if (!quiet) cli::cli_alert_success("Removed {length(removed)} diamondback temporary file{?s}.")
  invisible(removed)
}
