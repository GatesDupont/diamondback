#' Core area: habitat beyond a given distance from an edge
#'
#' Calculates, for each patch, the cells lying further than `depth` from the
#' patch edge, using an exact Euclidean distance transform
#' (`scipy.ndimage.distance_transform_edt`) rather than repeated erosion. One
#' pass, any depth, and non-square cells handled exactly via anisotropic
#' sampling.
#'
#' @section How distance is measured:
#' The transform gives, for each habitat cell, the Euclidean distance from its
#' **centre** to the **centre of the nearest edge-source cell**. So a cell
#' directly adjacent to non-habitat has a distance of one cell width, not zero:
#' the physical patch boundary lies half a cell away from either centre. Core
#' cells are those with distance **strictly greater** than `depth`.
#'
#' The practical consequence: with 30 m cells and `depth = 300`, a cell needs at
#' least 11 cells of habitat between it and the nearest edge cell to qualify. If
#' you want the distance measured to the physical boundary instead, pass
#' `depth + res/2`.
#'
#' @section What counts as an edge:
#' `edge = "all"` (default) treats every non-habitat cell as an edge source,
#' including masked-out cells, missing cells, and everything beyond the grid.
#' Patches clipped by the study-area boundary lose core area accordingly, which
#' is conservative and usually right when the boundary is a real edge such as a
#' coastline.
#'
#' `edge = "background"` treats only valid background as an edge source. Masked,
#' missing and off-grid cells count as habitat, so a patch running off the study
#' area is not penalised for it. This is the right choice when the study-area
#' boundary is an arbitrary analytical cut through continuous habitat.
#'
#' @section Geographic rasters:
#' The distance transform needs uniform cell spacing, which a lon/lat grid does
#' not have. This function errors on lon/lat unless `units = "cells"`. Project to
#' an equal-area CRS first.
#'
#' @param x A [patch_result], or a labelled `SpatRaster`.
#' @param depth Edge depth. In map units unless `units = "cells"`.
#' @param units Units of `depth` and of the returned core area columns: `"m"`
#'   (default), `"km"`, or `"cells"`.
#' @param edge What counts as an edge source: `"all"` or `"background"`.
#' @param core_raster Also return a raster of core cells, as
#'   `result$core_raster`, labelled with patch IDs and `0` elsewhere.
#' @param max_memory_frac,memory_limit Memory guards. The distance transform
#'   holds a float64 array, making this the most memory-hungry step in the
#'   package (roughly 13 bytes per cell).
#' @param quiet Suppress progress reporting.
#'
#' @return The [patch_result] with `core_cells`, `core_area_m2`, `core_area_ha`,
#'   `core_area_km2` (or `core_area_cells`) and `core_fraction` added to
#'   `$metrics`. `core_fraction` is core cells over total cells, so a patch with
#'   no interior gets `0`.
#'
#' @seealso [patch_metrics()]
#' @export
#' @examplesIf diamondback_ready()
#' # A 7x7 block of habitat: with 1 m cells and depth 1, only cells more than
#' # one cell from the edge survive, giving a 3x3 core.
#' m <- matrix(0, 9, 9)
#' m[2:8, 2:8] <- 1
#' res <- patch_core_area(label_patches(m, quiet = TRUE), depth = 1, quiet = TRUE)
#' res$metrics[, c("patch_id", "cells", "core_cells", "core_fraction")]
patch_core_area <- function(x,
                            depth,
                            units = c("m", "km", "cells"),
                            edge = c("all", "background"),
                            core_raster = FALSE,
                            max_memory_frac = 0.6,
                            memory_limit = NULL,
                            quiet = FALSE) {
  units <- match.arg(units)
  edge <- match.arg(edge)
  t0 <- Sys.time()

  if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) || depth < 0) {
    cli::cli_abort("{.arg depth} must be a single non-negative number.", call = NULL)
  }

  st <- db_result_state(x, quiet = quiet)
  r <- st$template
  n <- st$n
  lonlat <- db_is_lonlat(r)

  if (lonlat && !identical(units, "cells")) {
    cli::cli_abort(c(
      "Core area needs a projected raster.",
      "x" = "{.arg x} is in geographic (lon/lat) coordinates, where cell spacing varies with latitude.",
      "i" = "The Euclidean distance transform assumes uniform spacing, so a metric depth would be wrong.",
      "*" = "Project to an equal-area CRS, e.g. {.code terra::project(x, \"EPSG:5070\", method = \"near\")}.",
      "*" = 'Or measure in cells with {.code units = "cells"}, if that is meaningful for your grid.'
    ), call = NULL)
  }

  db_check_memory(terra::ncell(r), "core", max_memory_frac, memory_limit, quiet = quiet)

  # Sampling: physical size of one cell step in each axis. In "cells" units a
  # step is 1 by definition, which is what makes depth-in-cells work on lon/lat.
  sampling <- if (identical(units, "cells")) {
    c(1, 1)
  } else {
    c(terra::yres(r), terra::xres(r))  # (row step, col step)
  }
  depth_map <- if (identical(units, "km")) depth * 1000 else depth

  if (n == 0L) {
    st$result$metrics <- db_add_core_cols(st$result$metrics, numeric(0), numeric(0), units)
    return(st$result)
  }

  py <- db_py()
  classes <- st$classes
  n_classes <- if (is.null(classes)) 1L else length(classes)

  if (!quiet) {
    cli::cli_alert_info(
      "Computing core area (depth {depth} {units}, edge = {.val {edge}}) via distance transform ..."
    )
  }

  core_count <- numeric(n)
  core_mask <- NULL

  # One transform per class: a distance transform over "any patch" would treat
  # a neighbouring class as habitat, which is wrong for multi-class runs.
  pb <- db_progress("Distance transform", n_classes, quiet = quiet)
  for (k in seq_len(n_classes)) {
    out <- py_try(
      py$core_counts(st$code, st$labels, as.integer(n), as.integer(k - 1L),
                     depth_map, reticulate::tuple(sampling[1], sampling[2]), edge),
      "computing the distance transform"
    )
    core_count <- core_count + py_num(out, "core_count")[-1]
    if (isTRUE(core_raster)) {
      core_mask <- if (is.null(core_mask)) py_get(out, "core_mask") else {
        py_try(reticulate::import("numpy", convert = FALSE)$logical_or(
          core_mask, py_get(out, "core_mask")), "combining core masks")
      }
    }
    db_progress_step(pb, k)
  }
  db_progress_done(pb)

  met <- st$result$metrics
  if (is.null(met) || !nrow(met)) met <- data.frame(patch_id = seq_len(n))
  cells <- if ("cells" %in% names(met)) met$cells else {
    py_num(py_try(py$patch_stats(st$labels, as.integer(n)), "counting cells"), "count")[-1]
  }

  # Core area from cell count times cell area: every core cell is a whole cell,
  # so there is nothing to integrate here.
  cell_area <- if (identical(units, "cells")) 1 else terra::xres(r) * terra::yres(r)
  met <- db_add_core_cols(met, core_count, core_count * cell_area, units)
  met$core_fraction <- ifelse(cells > 0, core_count / cells, NA_real_)

  st$result$metrics <- met
  st$result$metadata$core <- list(depth = depth, units = units, edge = edge)

  if (isTRUE(core_raster)) {
    st$result$core_raster <- db_core_to_rast(py, st$labels, core_mask, r)
  }

  if (!quiet) {
    tot <- sum(core_count)
    cli::cli_alert_success(
      "Core area for {.val {n}} patch{?es} in {db_elapsed(t0)} \\
       ({.val {format(tot, big.mark = ',')}} core cell{?s})."
    )
  }
  st$result
}

db_add_core_cols <- function(met, core_count, core_area, units) {
  met$core_cells <- core_count
  if (identical(units, "cells")) {
    met$core_area_cells <- core_count
  } else {
    met$core_area_m2 <- core_area
    met$core_area_ha <- core_area / 1e4
    met$core_area_km2 <- core_area / 1e6
  }
  met
}

db_core_to_rast <- function(py, labels, core_mask, template) {
  np <- reticulate::import("numpy", convert = FALSE)
  core_labels <- py_try(np$where(core_mask, labels, np$asarray(0L)$astype(labels$dtype)),
                        "building the core raster")
  out <- terra::rast(template)
  nr <- terra::nrow(out); nc <- terra::ncol(out)
  vals <- integer(nr * nc)
  for (b in db_row_blocks(nr, nc)) {
    idx <- ((b$row - 1L) * nc + 1L):((b$row - 1L + b$nrows) * nc)
    rows <- py_try(np$asarray(core_labels[
      reticulate::import_builtins(convert = FALSE)$slice(
        as.integer(b$row - 1L), as.integer(b$row - 1L + b$nrows))
    ])$astype(np$int32)$ravel(), "extracting core cells")
    vals[idx] <- as.integer(reticulate::py_to_r(rows))
  }
  terra::values(out) <- vals
  names(out) <- "core_patch_id"
  out
}
