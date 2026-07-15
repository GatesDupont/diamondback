# Patch geometry. See DESIGN.md section 4 for the perimeter definition; the
# short version is that only shared cell edges count, corners contribute
# nothing, and boundary against the study area is reported separately from
# boundary against real habitat.

#' Calculate patch-level metrics
#'
#' Geometry for every patch in a labelled result: area, perimeter, edge
#' composition, bounding box and centroid. Separated from [label_patches()]
#' because labelling is the expensive, cacheable step while metrics are cheap
#' and often wanted several ways. [analyze_patches()] does both in one call.
#'
#' @section What perimeter means:
#' Perimeter is the total length of **shared cell edges** between a patch cell
#' and any cell that is not part of that patch, plus any edges on the grid
#' border. Diagonal contact contributes **nothing**: two cells touching only at
#' a corner share no physical boundary. This holds regardless of whether the
#' patch was labelled with 4- or 8-connectivity --- connectivity decides what is
#' one patch, not where its boundary runs.
#'
#' Non-square cells are handled correctly. An edge shared with a left/right
#' neighbour runs north--south and is `yres` long; an edge shared with an
#' up/down neighbour runs east--west and is `xres` long. The counts of each are
#' reported separately (`edge_ns_cells`, `edge_ew_cells`) as exact integers.
#'
#' @section Habitat edge versus domain edge:
#' `perimeter = edge_valid + edge_domain`, and the split matters:
#'
#' * `edge_valid` --- boundary against an in-domain, known cell. The patch
#'   genuinely stops here. This is real habitat edge.
#' * `edge_domain` --- boundary against a masked-out cell, a missing cell, or the
#'   grid border. The patch may well continue; we cannot see it.
#'
#' `touches_domain_edge` flags patches with any `edge_domain`. **Their `cells`,
#' `area_*` and `perimeter` are lower bounds**, and they are usually the patches
#' you want to exclude from a size distribution.
#'
#' @section Geographic rasters:
#' Lon/lat grids are supported exactly, not approximated as planar. Cell area
#' and east--west extent vary with latitude but are constant within a raster
#' row, so per-row geodesic geometry is used to weight the accumulation: `area_m2`
#' comes from [terra::cellSize()] and edge lengths from geodesic per-row
#' distances. Metres are the unit for lon/lat regardless of `units`.
#'
#' @param x A [patch_result] from [label_patches()], or a labelled `SpatRaster`.
#' @param units Map units for area and length columns: `"m"` (default), `"km"`,
#'   or `"cells"` for raw cell counts and cell-edge counts. Ignored for lon/lat
#'   rasters, which always report metres.
#' @param edge_depth Optional core-area depth. When given, core columns are added
#'   as if by [patch_core_area()], in `units`.
#' @param bbox Include bounding-box columns.
#' @param centroid Include centroid columns.
#' @param quiet Suppress progress reporting.
#'
#' @return The [patch_result] with `$metrics` replaced by a data frame with one
#'   row per patch and these columns:
#'
#'   * `patch_id`, `class` (multi-class runs only), `cells`
#'   * `area_m2`, `area_ha`, `area_km2` --- area in map units
#'   * `perimeter` --- total exposed boundary
#'   * `edge_valid`, `edge_domain` --- perimeter split by what is on the far side
#'   * `edge_ns_cells`, `edge_ew_cells` --- exact edge counts, CRS-free
#'   * `edge_area_ratio` --- `perimeter / area`
#'   * `touches_domain_edge` --- whether metrics are a lower bound
#'   * `xmin`, `xmax`, `ymin`, `ymax` --- bounding box in map coordinates
#'   * `x_centroid`, `y_centroid` --- mean cell centre, which for a concave or
#'     multi-lobed patch may fall outside the patch
#'   * `core_cells`, `core_area_*`, `core_fraction` --- when `edge_depth` is set
#'
#' @seealso [patch_core_area()], [analyze_patches()]
#' @export
#' @examples
#' m <- matrix(c(1, 1, 0, 0,
#'               1, 1, 0, 0,
#'               0, 0, 0, 1,
#'               0, 0, 0, 1), nrow = 4, byrow = TRUE)
#' res <- patch_metrics(label_patches(m, quiet = TRUE))
#' res$metrics[, c("patch_id", "cells", "perimeter", "touches_domain_edge")]
patch_metrics <- function(x,
                          units = c("m", "km", "cells"),
                          edge_depth = NULL,
                          bbox = TRUE,
                          centroid = TRUE,
                          quiet = FALSE) {
  units <- match.arg(units)
  t0 <- Sys.time()

  st <- db_result_state(x, quiet = quiet)
  py <- db_py()
  n <- st$n
  r <- st$template

  if (n == 0L) {
    st$result$metrics <- db_empty_metrics(units, !is.null(edge_depth), bbox, centroid,
                                          multiclass = !is.null(st$classes) && length(st$classes) > 1L)
    return(st$result)
  }

  lonlat <- db_is_lonlat(r)
  rg <- db_row_geometry(r)

  if (!quiet) cli::cli_alert_info("Computing patch statistics ...")
  stats <- py_try(
    py$patch_stats(st$labels, as.integer(n),
                   cell_area_by_row = if (lonlat) rg$cell_area_by_row else NULL),
    "computing patch statistics"
  )
  cells <- py_num(stats, "count")[-1]
  sum_row <- py_num(stats, "sum_row")[-1]
  sum_col <- py_num(stats, "sum_col")[-1]

  if (!quiet) cli::cli_alert_info("Measuring patch boundaries ...")
  edges <- py_try(
    py$edge_lengths(st$labels, st$code, as.integer(n),
                    rg$dy_by_row, rg$dx_by_gridline,
                    na_background = st$na_background),
    "measuring patch boundaries"
  )

  edge_valid <- py_num(edges, "edge_valid_len")[-1]
  edge_domain <- py_num(edges, "edge_domain_len")[-1]
  edge_ns <- py_num(edges, "edge_ns_cells")[-1]
  edge_ew <- py_num(edges, "edge_ew_cells")[-1]

  area_m2 <- if (lonlat) py_num(stats, "area")[-1] else cells * rg$cell_area

  out <- data.frame(patch_id = seq_len(n))
  if (!is.null(st$classes) && length(st$classes) > 1L) {
    out$class <- st$classes[st$patch_class]
  }
  out$cells <- cells

  # Unit handling. "cells" reports raw counts, which is the only honest option
  # for a bare matrix with no CRS; anything else needs real map units.
  if (identical(units, "cells")) {
    out$area_cells <- cells
    out$perimeter_cells <- edge_ns + edge_ew
    out$edge_valid_cells <- py_num(edges, "edge_valid_cells")[-1]
    out$edge_domain_cells <- py_num(edges, "edge_domain_cells")[-1]
    out$perimeter <- out$perimeter_cells
    out$edge_valid <- out$edge_valid_cells
    out$edge_domain <- out$edge_domain_cells
  } else {
    scale_len <- if (lonlat || identical(units, "m")) 1 else 1e-3
    scale_area <- if (lonlat || identical(units, "m")) 1 else 1e-6

    out$area_m2 <- area_m2
    out$area_ha <- area_m2 / 1e4
    out$area_km2 <- area_m2 / 1e6
    if (identical(units, "km")) {
      out$area <- area_m2 * scale_area
    }
    out$perimeter <- (edge_valid + edge_domain) * scale_len
    out$edge_valid <- edge_valid * scale_len
    out$edge_domain <- edge_domain * scale_len
  }

  out$edge_ns_cells <- edge_ns
  out$edge_ew_cells <- edge_ew
  out$edge_area_ratio <- out$perimeter / ifelse(out$cells > 0,
                                                if (identical(units, "cells")) out$area_cells else out$area_m2,
                                                NA_real_)
  out$touches_domain_edge <- edge_domain > 0

  if (isTRUE(bbox)) {
    if (!quiet) cli::cli_alert_info("Computing bounding boxes ...")
    bb <- py_try(py$patch_bboxes(st$labels, as.integer(n)), "computing bounding boxes")
    rmin <- py_num(bb, "row_min")[-1]; rmax <- py_num(bb, "row_max")[-1]
    cmin <- py_num(bb, "col_min")[-1]; cmax <- py_num(bb, "col_max")[-1]
    e <- terra::ext(r)
    out$xmin <- e$xmin + cmin * rg_xres(r)
    out$xmax <- e$xmin + (cmax + 1) * rg_xres(r)
    # Row 0 is the top row, so ymax comes from row_min.
    out$ymax <- e$ymax - rmin * rg_yres(r)
    out$ymin <- e$ymax - (rmax + 1) * rg_yres(r)
  }

  if (isTRUE(centroid)) {
    e <- terra::ext(r)
    mean_col <- sum_col / cells
    mean_row <- sum_row / cells
    out$x_centroid <- e$xmin + (mean_col + 0.5) * rg_xres(r)
    out$y_centroid <- e$ymax - (mean_row + 0.5) * rg_yres(r)
  }

  st$result$metrics <- out
  st$result$metadata$metrics_units <- units

  if (!is.null(edge_depth)) {
    st$result <- patch_core_area(st$result, depth = edge_depth, units = units, quiet = quiet)
  }

  if (!quiet) cli::cli_alert_success("Metrics for {.val {n}} patch{?es} in {db_elapsed(t0)}.")
  st$result
}

rg_xres <- function(r) terra::xres(r)
rg_yres <- function(r) terra::yres(r)

db_empty_metrics <- function(units, core, bbox, centroid, multiclass) {
  base <- data.frame(patch_id = integer(), cells = numeric())
  if (multiclass) base$class <- numeric()
  cols <- if (identical(units, "cells")) {
    c("area_cells", "perimeter", "edge_valid", "edge_domain")
  } else {
    c("area_m2", "area_ha", "area_km2", "perimeter", "edge_valid", "edge_domain")
  }
  cols <- c(cols, "edge_ns_cells", "edge_ew_cells", "edge_area_ratio")
  for (cn in cols) base[[cn]] <- numeric()
  base$touches_domain_edge <- logical()
  if (bbox) for (cn in c("xmin", "xmax", "ymin", "ymax")) base[[cn]] <- numeric()
  if (centroid) for (cn in c("x_centroid", "y_centroid")) base[[cn]] <- numeric()
  if (core) for (cn in c("core_cells", "core_fraction")) base[[cn]] <- numeric()
  base
}

#' Resolve a result (or bare raster) into the arrays the kernels need
#'
#' Reuses the arrays from the labelling run when they are still alive, and
#' rebuilds them from the labelled raster otherwise --- which is what makes a
#' result read back from disk work exactly like a fresh one.
#' @noRd
db_result_state <- function(x, quiet = FALSE) {
  if (inherits(x, "patch_result")) {
    arrays <- .subset2(x, "arrays")
    if (!is.null(arrays) && db_array_alive(arrays$labels) && db_array_alive(arrays$code)) {
      return(list(
        result = x, labels = arrays$labels, code = arrays$code, n = arrays$n,
        classes = arrays$classes, na_background = isTRUE(arrays$na_background),
        patch_class = if (!is.null(x$metrics$class)) {
          match(x$metrics$class, arrays$classes)
        } else NULL,
        template = db_patches(x)
      ))
    }
    if (!quiet) cli::cli_alert_info("Rebuilding arrays from the labelled raster ...")
    return(db_state_from_labels(db_patches(x), result = x, quiet = quiet))
  }

  if (inherits(x, "SpatRaster") || is.matrix(x) || is.character(x)) {
    r <- db_as_rast(x)
    return(db_state_from_labels(r, result = NULL, quiet = quiet))
  }

  cli::cli_abort(c(
    "{.arg x} must be a {.cls patch_result} or a labelled {.cls SpatRaster}.",
    "x" = "Got {.cls {class(x)[1]}}."
  ), call = NULL)
}

#' Rebuild code and label arrays from a labelled raster
#'
#' A labelled raster carries three of the four states (patch / background / NA)
#' but cannot distinguish outside-domain from missing, so NA is read back as
#' missing. That is the conservative choice: missing and outside behave
#' identically for edge classification, so boundary against either is still
#' counted as domain edge rather than habitat edge.
#' @noRd
db_state_from_labels <- function(r, result = NULL, quiet = FALSE) {
  py <- db_py()
  np <- reticulate::import("numpy", convert = FALSE)
  nr <- terra::nrow(r); nc <- terra::ncol(r)

  db_check_memory(terra::ncell(r), "metrics", quiet = TRUE)

  use_int64 <- terra::ncell(r) > .Machine$integer.max
  labels <- py_try(np$zeros(reticulate::tuple(as.integer(nr), as.integer(nc)),
                            dtype = if (use_int64) np$int64 else np$int32),
                   "allocating the label array")
  code <- py_try(np$zeros(reticulate::tuple(as.integer(nr), as.integer(nc)), dtype = np$uint8),
                 "allocating the cell-state array")

  blocks <- db_row_blocks(nr, nc)
  terra::readStart(r); on.exit(terra::readStop(r), add = TRUE)
  pb <- db_progress("Reading labels", length(blocks), quiet = quiet)
  n_max <- 0
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    v <- terra::readValues(r, row = b$row, nrows = b$nrows, col = 1, ncols = nc, mat = FALSE)
    v <- as.numeric(v)
    n_max <- max(n_max, suppressWarnings(max(v, na.rm = TRUE)), na.rm = TRUE)
    # Label raster convention: NA -> outside, 0 -> background, >0 -> foreground.
    lab_blk <- py_try(py$code_block(v, as.integer(b$nrows), as.integer(nc),
                                    mask = NULL, class_values = NULL),
                      "converting labels to cell states")
    py_try(py$set_rows(code, as.integer(b$row - 1L),
                       as.integer(b$row - 1L + b$nrows), lab_blk),
           "filling the cell-state array")
    lv <- v; lv[is.na(lv)] <- 0
    blk <- py_try(np$asarray(matrix(lv, nrow = b$nrows, ncol = nc, byrow = TRUE))$astype(labels$dtype),
                  "converting label values")
    py_try(py$set_rows(labels, as.integer(b$row - 1L),
                       as.integer(b$row - 1L + b$nrows), blk),
           "filling the label array")
    db_progress_step(pb, i)
  }
  db_progress_done(pb)

  if (!is.finite(n_max)) n_max <- 0
  n <- as.integer(n_max)

  res <- result
  if (is.null(res)) {
    counts <- db_code_counts(code)
    res <- new_patch_result(
      patches = if (terra::inMemory(r)) r else NULL,
      patches_path = if (terra::inMemory(r)) NULL else terra::sources(r)[1],
      metrics = data.frame(patch_id = seq_len(n), cells = rep(NA_real_, n)),
      metadata = db_metadata(
        source = db_source_info(r), mask_source = NULL, geometry = db_geometry(r),
        class = NULL, directions = NA_integer_, na = "outside", crop = FALSE,
        counts = counts, n_patches = n, elapsed = NA_real_
      )
    )
  }

  list(result = res, labels = labels, code = code, n = n,
       classes = NULL, na_background = FALSE, patch_class = NULL,
       template = r)
}
