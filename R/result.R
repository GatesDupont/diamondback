# The patch_result S3 class: a list with an active $patches accessor, so
# callers cannot tell whether the raster is in memory or on disk.

#' Patch analysis results
#'
#' The object returned by [label_patches()] and [analyze_patches()]. It is a
#' plain list underneath, so nothing is hidden, with three components you will
#' use directly:
#'
#' * `$patches` --- the labelled `SpatRaster`. This is an accessor, not stored
#'   data: if the run wrote to disk it opens the file; if it stayed in memory it
#'   returns the raster it holds. Either way you get a `SpatRaster`.
#' * `$metrics` --- one row per patch. After [label_patches()] this has patch ID
#'   and cell count; after [patch_metrics()] or [analyze_patches()] it has the
#'   full geometry.
#' * `$metadata` --- everything needed to reproduce the run: source fingerprint,
#'   geometry, class, connectivity, NA policy, cell-state counts, and backend
#'   versions.
#'
#' `$path` gives the file the labels were written to, or `NULL`.
#'
#' @section Arrays:
#' A result also carries a handle to the NumPy arrays from the run, which is
#' what lets [patch_metrics()] and [compare_patches()] avoid relabelling. These
#' do not survive [saveRDS()] or a session restart; the arrays are rebuilt from
#' the labelled raster when needed, so a restored result still works.
#'
#' @name patch_result
#' @seealso [label_patches()], [patch_metrics()], [write_patch_result()]
NULL

new_patch_result <- function(patches, patches_path, metrics, metadata, arrays = NULL) {
  structure(
    list(
      patches_obj = patches,
      path = patches_path,
      metrics = metrics,
      metadata = metadata,
      arrays = arrays
    ),
    class = "patch_result"
  )
}

#' @export
`$.patch_result` <- function(x, name) {
  if (identical(name, "patches")) {
    return(db_patches(x))
  }
  NextMethod()
}

#' @export
`[[.patch_result` <- function(x, i, ...) {
  if (identical(i, "patches")) {
    return(db_patches(x))
  }
  NextMethod()
}

#' Get the labelled raster from a result
#' @noRd
db_patches <- function(x) {
  obj <- .subset2(x, "patches_obj")
  if (!is.null(obj)) {
    # A SpatRaster whose C++ pointer died with the session is unusable; catch
    # that here rather than letting terra fail confusingly later.
    ok <- tryCatch({ terra::ncell(obj); TRUE }, error = function(e) FALSE)
    if (ok) return(obj)
  }
  p <- .subset2(x, "path")
  if (!is.null(p) && file.exists(p)) return(terra::rast(p))
  if (!is.null(p)) {
    cli::cli_abort(c(
      "The labelled raster for this result is missing.",
      "x" = "Expected it at {.path {p}}.",
      "i" = "It may have been moved or deleted since the analysis ran."
    ), call = NULL)
  }
  cli::cli_abort(c(
    "This result no longer holds a labelled raster.",
    "i" = "In-memory rasters do not survive {.fn saveRDS} or a session restart.",
    "i" = "Use {.fn write_patch_result} to save results durably."
  ), call = NULL)
}

#' @export
print.patch_result <- function(x, ...) {
  m <- .subset2(x, "metadata")
  g <- m$geometry
  cli::cli_h1("<patch_result>")

  # qty() is needed because the interpolated value is a formatted string; cli
  # would otherwise infer the quantity from it and print "0 patch".
  cli::cli_text("{.strong {format(m$n_patches, big.mark = ',')}} \\
                 {cli::qty(m$n_patches)}patch{?es} \\
                 from {.val {format(g$ncell, big.mark = ',')}} cells \\
                 ({g$nrow} x {g$ncol})")

  src <- if (identical(m$source$type, "file")) basename(m$source$path) else paste0("<", m$source$type, ">")
  cli::cli_dl(c(
    "source" = src,
    "class" = if (isTRUE(m$binary)) "binary (non-zero = foreground)" else paste(m$class, collapse = ", "),
    "connectivity" = paste0(m$directions, "-neighbour"),
    "NA policy" = m$na,
    "resolution" = sprintf("%g x %g%s", g$xres, g$yres, if (isTRUE(g$lonlat)) " (lon/lat)" else "")
  ))

  cl <- m$cells
  cli::cli_text("")
  cli::cli_text("{.strong Cells}")
  cli::cli_dl(c(
    "patch" = format(cl$foreground, big.mark = ","),
    "background" = format(cl$background, big.mark = ","),
    "missing" = format(cl$missing, big.mark = ","),
    "outside domain" = format(cl$outside, big.mark = ",")
  ))

  met <- .subset2(x, "metrics")
  cli::cli_text("")
  has_geom <- "area_m2" %in% names(met)
  cli::cli_text("{.strong $metrics}: {nrow(met)} row{?s} x {ncol(met)} column{?s}\\
                 {if (!has_geom) ' (cell counts only; run patch_metrics() for geometry)' else ''}")

  p <- .subset2(x, "path")
  if (!is.null(p)) {
    cli::cli_text("{.strong $patches}: {.path {p}}")
  } else {
    cli::cli_text("{.strong $patches}: in memory")
  }
  if (isTRUE(m$validated)) cli::cli_alert_success("Validated.")
  invisible(x)
}

#' @export
summary.patch_result <- function(object, ...) {
  m <- .subset2(object, "metadata")
  met <- .subset2(object, "metrics")
  print(object)

  if (nrow(met) == 0) {
    cli::cli_alert_warning("No patches to summarise.")
    return(invisible(list(metadata = m, metrics = met)))
  }

  cli::cli_h2("Patch size")
  q <- stats::quantile(met$cells, c(0, 0.5, 0.9, 1))
  cli::cli_dl(c(
    "cells (min / median / p90 / max)" =
      paste(format(round(q), big.mark = ","), collapse = " / "),
    "total patch cells" = format(sum(met$cells), big.mark = ",")
  ))

  if ("area_ha" %in% names(met)) {
    qa <- stats::quantile(met$area_ha, c(0, 0.5, 0.9, 1))
    cli::cli_dl(c(
      "area ha (min / median / p90 / max)" =
        paste(signif(qa, 4), collapse = " / "),
      "total area (ha)" = signif(sum(met$area_ha), 6)
    ))
  }
  if ("touches_domain_edge" %in% names(met)) {
    n_edge <- sum(met$touches_domain_edge)
    cli::cli_text("")
    cli::cli_alert_info(
      "{n_edge} patch{?es} ({round(100 * n_edge / nrow(met), 1)}%) touch the domain edge; \\
       their area and edge metrics are lower bounds."
    )
  }
  if ("class" %in% names(met)) {
    cli::cli_h2("Patches per class")
    tb <- table(met$class)
    cli::cli_dl(stats::setNames(as.character(as.integer(tb)), names(tb)))
  }
  invisible(list(metadata = m, metrics = met))
}

#' Reconstruct the four-state cell domain as a raster
#'
#' The labelled raster collapses "outside the domain" and "missing" into `NA`,
#' because a raster has only one no-data value. This rebuilds the full
#' distinction as a categorical raster with levels `outside`, `missing`,
#' `background` and `patch`.
#'
#' The source raster is re-read, so this costs a pass over the data. The exact
#' counts are always available without that, in `result$metadata$cells`.
#'
#' @param x A [patch_result].
#' @param filename Optional path to write to.
#' @param overwrite Overwrite `filename` if it exists.
#' @return A categorical `SpatRaster`.
#' @export
#' @examplesIf diamondback_ready()
#' m <- matrix(c(1, 1, NA, 0, 0, 1, 0, 1, 1), nrow = 3)
#' res <- label_patches(m, quiet = TRUE)
#' terra::freq(patch_domain(res))
patch_domain <- function(x, filename = NULL, overwrite = FALSE) {
  stopifnot(inherits(x, "patch_result"))
  p <- db_patches(x)
  arrays <- .subset2(x, "arrays")

  code <- if (!is.null(arrays) && db_array_alive(arrays$code)) {
    arrays$code
  } else {
    cli::cli_abort(c(
      "The cell-state array for this result is no longer available.",
      "i" = "It does not survive a session restart or {.fn saveRDS}.",
      "i" = "Re-run {.fn label_patches} to rebuild it."
    ), call = NULL)
  }

  py <- db_py()
  out <- terra::rast(p)
  nr <- terra::nrow(out); nc <- terra::ncol(out)
  vals <- integer(nr * nc)
  for (b in db_row_blocks(nr, nc)) {
    rows <- py_try(py$domain_rows(code, as.integer(b$row - 1L),
                                  as.integer(b$row - 1L + b$nrows)),
                   "extracting cell states")
    idx <- ((b$row - 1L) * nc + 1L):((b$row - 1L + b$nrows) * nc)
    vals[idx] <- as.integer(py_r(rows))
  }
  # Every foreground class collapses to a single "patch" level: the point of
  # this raster is the four states, and class identity is in $metrics.
  vals[vals >= 3L] <- 3L
  terra::values(out) <- vals
  levels(out) <- data.frame(
    value = 0:3,
    domain = c("outside", "missing", "background", "patch")
  )
  names(out) <- "domain"
  if (!is.null(filename)) {
    db_register_file(filename)
    out <- terra::writeRaster(out, filename, overwrite = overwrite)
  }
  out
}

# A NumPy array handle from a dead session errors on access; test cheaply.
db_array_alive <- function(a) {
  if (is.null(a)) return(FALSE)
  tryCatch({ reticulate::py_to_r(a$shape); TRUE }, error = function(e) FALSE)
}
