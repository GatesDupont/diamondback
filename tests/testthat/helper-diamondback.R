# Every test needs the Python backend. Rather than let each test fail
# separately on a machine without it, skip once with a clear reason.
skip_if_no_python <- function() {
  ok <- tryCatch({
    diamondback_python(quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(ok, "Python backend (NumPy + SciPy) not available")
}

# A projected raster with 10 m square cells, built from a matrix so that tests
# can state the expected answer as a picture.
rast_from <- function(m, res = 10, crs = "EPSG:5070", yres = res) {
  r <- terra::rast(
    nrows = nrow(m), ncols = ncol(m),
    xmin = 0, xmax = ncol(m) * res,
    ymin = 0, ymax = nrow(m) * yres,
    crs = crs
  )
  terra::values(r) <- as.vector(t(m))
  r
}

# terra::patches() wants NA background; diamondback wants 0. This bridges them
# so the two can be compared on the same input.
terra_patch_count <- function(m, directions = 8) {
  r <- rast_from(m)
  r[r == 0] <- NA
  p <- terra::patches(r, directions = directions, zeroAsNA = FALSE)
  length(unique(stats::na.omit(terra::values(p, mat = FALSE))))
}

# Two labellings agree if they induce the same partition of cells, regardless
# of the order the labels were handed out in.
same_partition <- function(a, b) {
  a <- as.vector(a); b <- as.vector(b)
  a[is.na(a)] <- 0; b[is.na(b)] <- 0
  fa <- a > 0; fb <- b > 0
  if (!identical(fa, fb)) return(FALSE)
  if (!any(fa)) return(TRUE)
  # A bijection between label sets exists iff the cross-tabulation has exactly
  # one non-zero entry per row and per column.
  tb <- table(a[fa], b[fb])
  all(rowSums(tb > 0) == 1) && all(colSums(tb > 0) == 1)
}
