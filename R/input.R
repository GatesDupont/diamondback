# Normalising inputs into a "grid": everything downstream needs to know about
# the geometry of an analysis, resolved once and validated once.

#' Coerce raster-like input to a single-layer SpatRaster
#'
#' Accepts a SpatRaster, a filename, or a matrix. A matrix is given a trivial
#' 1-unit-per-cell extent and no CRS, which is what the tests and small examples
#' want; anything spatial should come in as a raster.
#'
#' @noRd
db_as_rast <- function(x, arg = "x") {
  if (inherits(x, "SpatRaster")) {
    r <- x
  } else if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      cli::cli_abort(c(
        "{.arg {arg}} looks like a filename but the file does not exist.",
        "x" = "{.path {x}}"
      ), call = NULL)
    }
    r <- terra::rast(x)
  } else if (is.matrix(x)) {
    r <- terra::rast(
      nrows = nrow(x), ncols = ncol(x),
      xmin = 0, xmax = ncol(x), ymin = 0, ymax = nrow(x),
      crs = ""
    )
    terra::values(r) <- as.vector(t(x))
  } else if (inherits(x, "patch_result")) {
    cli::cli_abort(c(
      "{.arg {arg}} is a {.cls patch_result}; this function expects a raster.",
      "i" = "Use {.code {arg}$patches} for the labelled raster."
    ), call = NULL)
  } else {
    cli::cli_abort(c(
      "{.arg {arg}} must be a {.cls SpatRaster}, a raster filename, or a matrix.",
      "x" = "Got {.cls {class(x)[1]}}."
    ), call = NULL)
  }

  if (terra::nlyr(r) > 1L) {
    cli::cli_abort(c(
      "{.arg {arg}} has {terra::nlyr(r)} layers; diamondback works on one layer at a time.",
      "i" = "Select a layer, e.g. {.code {arg}[[1]]}, or use {.fn track_patch_series} for a time series."
    ), call = NULL)
  }
  if (terra::ncell(r) == 0) {
    cli::cli_abort("{.arg {arg}} has no cells.", call = NULL)
  }
  r
}

#' Is this raster in geographic coordinates?
#'
#' Deliberately does not use `perhaps = TRUE`, which guesses from the extent: a
#' CRS-less 4x4 matrix looks like degrees under that rule and silently acquires
#' geodesic perimeters in the hundreds of kilometres. An unknown CRS is treated
#' as planar, and its map units are taken at face value.
#' @noRd
db_is_lonlat <- function(r) {
  # terra warns about an unknown CRS even with warn = FALSE; an unknown CRS is
  # an expected, handled case here, so the warning is noise.
  isTRUE(suppressWarnings(terra::is.lonlat(r, warn = FALSE)))
}

#' Describe the geometry of a raster in a comparable, serialisable form
#' @noRd
db_geometry <- function(r) {
  e <- terra::ext(r)
  # Unnamed doubles throughout. terra::ext() returns *named* numerics, and a
  # JSON round-trip turns a whole number into an integer -- either one makes an
  # otherwise identical geometry compare unequal, which would silently defeat
  # the cache. Normalising here means the comparison stays a plain all.equal().
  num <- function(v) as.numeric(unname(v))
  list(
    nrow = num(terra::nrow(r)),
    ncol = num(terra::ncol(r)),
    ncell = num(terra::ncell(r)),
    xmin = num(e$xmin), xmax = num(e$xmax),
    ymin = num(e$ymin), ymax = num(e$ymax),
    xres = num(terra::xres(r)), yres = num(terra::yres(r)),
    crs = terra::crs(r),
    lonlat = db_is_lonlat(r)
  )
}

#' Compare two geometries, returning the names of fields that differ
#' @noRd
db_geometry_diff <- function(g1, g2, tol = 1e-6) {
  diffs <- character()
  num_eq <- function(a, b) isTRUE(abs(a - b) <= tol * max(1, abs(a), abs(b)))

  if (!identical(g1$nrow, g2$nrow) || !identical(g1$ncol, g2$ncol)) {
    diffs <- c(diffs, "dimensions")
  }
  if (!num_eq(g1$xres, g2$xres) || !num_eq(g1$yres, g2$yres)) {
    diffs <- c(diffs, "resolution")
  }
  if (!num_eq(g1$xmin, g2$xmin) || !num_eq(g1$xmax, g2$xmax) ||
      !num_eq(g1$ymin, g2$ymin) || !num_eq(g1$ymax, g2$ymax)) {
    diffs <- c(diffs, "extent")
  }
  # Origin is the offset of the grid from a multiple of the resolution; two
  # rasters can share resolution and differ here, which silently misaligns cells.
  o1 <- c(g1$xmin %% g1$xres, g1$ymin %% g1$yres)
  o2 <- c(g2$xmin %% g2$xres, g2$ymin %% g2$yres)
  if (!num_eq(o1[1], o2[1]) || !num_eq(o1[2], o2[2])) {
    diffs <- c(diffs, "origin")
  }
  if (!db_crs_same(g1$crs, g2$crs)) {
    diffs <- c(diffs, "CRS")
  }
  unique(diffs)
}

db_crs_same <- function(a, b) {
  if (identical(a, b)) return(TRUE)
  # Empty CRS on both sides (bare matrices) counts as compatible.
  if (!nzchar(trimws(a)) && !nzchar(trimws(b))) return(TRUE)
  ok <- tryCatch(terra::same.crs(a, b), error = function(e) NA)
  isTRUE(ok)
}

#' Per-row geometry vectors used to weight edge and area accumulation
#'
#' For a projected raster these are constant. For lon/lat they vary with
#' latitude but are constant within a row, which is what makes the weighted
#' bincount approach exact rather than an approximation of a variable field.
#'
#' Returns:
#'   dy_by_row     length nrow, north-south extent of a cell in each row
#'   dx_by_gridline length nrow+1, east-west extent at each horizontal grid
#'                  line (element j is the line above row j)
#'   cell_area_by_row length nrow, or NULL when area is a constant
#' @noRd
db_row_geometry <- function(r) {
  nr <- terra::nrow(r)
  if (!db_is_lonlat(r)) {
    return(list(
      dy_by_row = rep(terra::yres(r), nr),
      dx_by_gridline = rep(terra::xres(r), nr + 1L),
      cell_area_by_row = NULL,
      cell_area = terra::xres(r) * terra::yres(r)
    ))
  }

  e <- terra::ext(r)
  xres <- terra::xres(r)
  yres <- terra::yres(r)

  # Latitudes of the nrow+1 horizontal grid lines, top to bottom.
  lat_lines <- seq(e$ymax, e$ymin, length.out = nr + 1L)
  lat_centres <- (lat_lines[-1] + lat_lines[-(nr + 1L)]) / 2

  # East-west extent of one cell at a given latitude, geodesic.
  dx_at <- function(lat) {
    p1 <- cbind(e$xmin, lat)
    p2 <- cbind(e$xmin + xres, lat)
    terra::distance(p1, p2, lonlat = TRUE, pairwise = TRUE)
  }
  dx_by_gridline <- as.numeric(dx_at(lat_lines))

  # North-south extent of a cell: geodesic and near-constant, but computed per
  # row rather than assumed.
  dy_by_row <- as.numeric(terra::distance(
    cbind(e$xmin, lat_lines[-(nr + 1L)]),
    cbind(e$xmin, lat_lines[-1]),
    lonlat = TRUE, pairwise = TRUE
  ))

  # Exact geodesic cell area; constant within a row, so one column suffices.
  cs <- terra::cellSize(r, unit = "m", transform = FALSE)
  area_col <- terra::values(cs[, 1, drop = FALSE], mat = FALSE)
  cell_area_by_row <- as.numeric(area_col)[seq_len(nr)]

  list(
    dy_by_row = dy_by_row,
    dx_by_gridline = dx_by_gridline,
    cell_area_by_row = cell_area_by_row,
    cell_area = NA_real_
  )
}

#' Validate and normalise the `class` argument
#' @noRd
db_check_class <- function(class, x) {
  if (is.null(class)) return(NULL)
  if (identical(class, "all")) return("all")
  # Check for NA before checking the type: a bare `class = NA` is logical, and
  # "must be numeric" would be a much less useful thing to tell that user than
  # "NA is not a class".
  if (anyNA(class)) {
    cli::cli_abort(c(
      "{.arg class} contains {.val NA}.",
      "i" = "{.val NA} marks missing data in diamondback and can never be a patch class.",
      "i" = 'To treat {.val NA} as background, use {.code na = "background"}.'
    ), call = NULL)
  }
  if (!is.numeric(class)) {
    cli::cli_abort(c(
      '{.arg class} must be numeric, {.val NULL} (binary), or {.val "all"}.',
      "x" = "Got {.cls {class(class)[1]}}."
    ), call = NULL)
  }
  if (anyDuplicated(class)) {
    cli::cli_abort("{.arg class} contains duplicate values.", call = NULL)
  }
  as.numeric(class)
}

#' Resolve `class = "all"` into the actual values present
#' @noRd
db_discover_classes <- function(r, mask = NULL, quiet = FALSE) {
  if (!quiet) cli::cli_alert_info("Scanning for distinct raster values ({.arg class = \"all\"}) ...")
  u <- terra::unique(r, incomparables = TRUE)
  vals <- sort(unique(stats::na.omit(as.numeric(u[[1]]))))
  if (length(vals) == 0) {
    cli::cli_abort("{.arg x} contains no non-NA values, so there is nothing to label.", call = NULL)
  }
  if (length(vals) > 252) {
    cli::cli_abort(c(
      '{.arg class = "all"} found {length(vals)} distinct values; at most 252 are supported.',
      "i" = "This raster looks continuous rather than categorical.",
      "i" = "Pass the specific values you want via {.arg class}."
    ), call = NULL)
  }
  vals
}

#' Validate a mask raster against the target grid
#' @noRd
db_check_mask <- function(mask, r) {
  if (is.null(mask)) return(NULL)
  m <- db_as_rast(mask, arg = "mask")
  d <- db_geometry_diff(db_geometry(m), db_geometry(r))
  if (length(d)) {
    cli::cli_abort(c(
      "{.arg mask} does not align with {.arg x}.",
      "x" = "Mismatched: {.val {d}}.",
      "i" = "diamondback will not resample silently. Align the mask yourself, e.g. \\
             {.code mask <- terra::resample(mask, x, method = \"near\")}."
    ), call = NULL)
  }
  m
}
