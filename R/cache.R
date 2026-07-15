# Caching. The rule is that a cache hit must be provably the same computation,
# and a miss must say why. See DESIGN.md section 5.

#' Fields that define the identity of a labelling run
#'
#' Deliberately excludes output paths, progress settings and metric options:
#' none of them can change the labels, so none of them should throw away an
#' hour of work.
#'
#' Every leaf is canonicalised to a character scalar. The reason is that one
#' side of a cache comparison has been through JSON and the other has not, and
#' JSON does not preserve R's type distinctions: `0.0` comes back as an integer
#' and `NA_real_` does not come back as `NA_real_`. Comparing strings that both
#' sides generate the same way removes an entire family of false cache misses,
#' which are worse than they sound --- a cache that never hits is a cache that
#' silently costs an hour per run.
#' @noRd
db_cache_fields <- function(meta) {
  s <- function(v) {
    if (is.null(v) || length(v) == 0) return("<none>")
    if (all(is.na(v))) return("<na>")
    paste(format(v, digits = 15, trim = TRUE), collapse = "|")
  }
  src <- function(x) {
    if (is.null(x)) return(list(type = "<none>"))
    list(type = s(x$type), path = s(x$path), size = s(x$size),
         mtime = s(x$mtime), hash = s(x$hash))
  }
  g <- meta$geometry

  list(
    algorithm_version = s(meta$algorithm_version),
    source = src(meta$source),
    mask_source = src(meta$mask_source),
    geometry = list(
      nrow = s(g$nrow), ncol = s(g$ncol),
      xmin = s(g$xmin), xmax = s(g$xmax), ymin = s(g$ymin), ymax = s(g$ymax),
      xres = s(g$xres), yres = s(g$yres), crs = s(g$crs)
    ),
    # "binary" and "class" are one fact, so they are one key: NA-as-class does
    # not survive JSON, but the string "<binary>" does.
    class = if (isTRUE(meta$binary) || is.null(meta$class) || all(is.na(meta$class))) {
      "<binary>"
    } else {
      s(sort(as.numeric(meta$class)))
    },
    directions = s(meta$directions),
    na = s(meta$na),
    crop = s(meta$crop)
  )
}

db_cache_key <- function(meta) {
  digest::digest(db_cache_fields(meta), algo = "sha1")
}

#' Decide whether a cached result matches the request
#'
#' Returns TRUE, or FALSE with the mismatching field attached, so the caller can
#' say *why* it is redoing the work rather than silently redoing it.
#' @noRd
db_cache_match <- function(cached_meta, request_meta) {
  # Prefer the canonical fields stored at write time. They are strings, so JSON
  # returns them byte-identical; re-deriving them from the *parsed* metadata
  # would go through NA and numeric round-trips that JSON does not preserve.
  # The re-derivation is kept only for results written by older versions.
  a <- cached_meta$cache_fields %||% db_cache_fields(cached_meta)
  b <- db_cache_fields(request_meta)

  # Both sides are now canonical strings, so the comparison is exact equality
  # and the walk exists only to name the field that differs.
  cmp <- function(x, y, path = "") {
    if (is.list(x) || is.list(y)) {
      keys <- union(names(x), names(y))
      for (k in keys) {
        r <- cmp(x[[k]], y[[k]], paste0(path, if (nzchar(path)) "$", k))
        if (!isTRUE(r)) return(r)
      }
      return(TRUE)
    }
    if (identical(as.character(x), as.character(y))) return(TRUE)
    path
  }

  # An unfingerprintable source (a large in-memory raster) can never match: we
  # have no way to know it did not change. Better a slow rerun than a wrong one.
  if (identical(a$source$type, "memory") && identical(a$source$hash, "<na>")) {
    return("source (in-memory raster, cannot be fingerprinted)")
  }

  cmp(a, b)
}

#' Look for a reusable result in a directory
#' @noRd
db_cache_lookup <- function(output_dir, request_meta, quiet = FALSE) {
  f_meta <- file.path(output_dir, RESULT_FILES$metadata)
  f_lab <- file.path(output_dir, RESULT_FILES$labels)
  if (!file.exists(f_meta) || !file.exists(f_lab)) return(NULL)

  cached <- tryCatch(jsonlite::fromJSON(f_meta, simplifyVector = TRUE),
                     error = function(e) NULL)
  if (is.null(cached)) {
    if (!quiet) cli::cli_alert_warning("Cached metadata in {.path {output_dir}} is unreadable; recomputing.")
    return(NULL)
  }

  m <- db_cache_match(cached, request_meta)
  if (isTRUE(m)) {
    res <- tryCatch(read_patch_result(output_dir), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    if (!quiet) {
      cli::cli_alert_success(
        "Reusing cached result from {.path {output_dir}} \\
         ({format(cached$n_patches, big.mark = ',')} patches, labelled {cached$created})."
      )
    }
    return(res)
  }

  if (!quiet) {
    cli::cli_alert_info(
      "Cached result in {.path {output_dir}} does not match this request \\
       ({.field {m}} differs); recomputing."
    )
  }
  NULL
}

#' Build the metadata for a request without doing the work
#'
#' The cache has to be checked *before* labelling, which means the identifying
#' metadata has to be constructible from the inputs alone.
#' @noRd
db_request_meta <- function(x, class, directions, mask, na, crop) {
  r <- db_as_rast(x)
  mask_r <- db_check_mask(mask, r)
  if (is.null(crop)) crop <- !is.null(mask_r)
  if (isTRUE(crop) && !is.null(mask_r)) {
    r <- db_crop_to_mask(r, mask_r, quiet = TRUE)$x
  }
  cls <- db_check_class(class, r)
  if (identical(cls, "all")) cls <- db_discover_classes(r, quiet = TRUE)

  list(
    algorithm_version = ALGORITHM_VERSION,
    source = db_source_info(x),
    mask_source = if (is.null(mask)) NULL else db_source_info(mask),
    geometry = db_geometry(r),
    class = if (is.null(cls)) NA_real_ else cls,
    binary = is.null(cls),
    directions = directions,
    na = na,
    crop = crop
  )
}
