# The tests assume the backend already exists; they never provision it. On CI it
# is installed by the workflow (see .github/workflows/R-CMD-check.yaml), which
# pins the versions and puts them in the log, and RETICULATE_PYTHON points here.
# Locally, any Python with NumPy and SciPy will do, or run
# diamondback_install_python() once.
#
# Nothing here sets options(diamondback.python = "managed"). A test helper that
# authorises downloads is a config people copy, and it hides an install inside
# what looks like a test run. diamondback_install_python() is covered by its own
# CI job instead.
skip_if_no_python <- function() {
  err <- NULL
  ok <- tryCatch({
    diamondback_python(quiet = TRUE)
    TRUE
  }, error = function(e) {
    err <<- conditionMessage(e)
    FALSE
  })

  # Skipping is the right answer on a contributor's laptop and the wrong one on
  # CI, where every test skipping would still be a green tick. When
  # DIAMONDBACK_REQUIRE_BACKEND is set the environment is promising a backend,
  # so its absence is a failure rather than a shrug.
  if (!ok && nzchar(Sys.getenv("DIAMONDBACK_REQUIRE_BACKEND"))) {
    stop("DIAMONDBACK_REQUIRE_BACKEND is set, so the Python backend was expected ",
         "to be available, but it is not. RETICULATE_PYTHON=",
         encodeString(Sys.getenv("RETICULATE_PYTHON"), quote = '"'),
         ". Original error: ", err %||% "(none)")
  }

  testthat::skip_if_not(
    ok,
    paste("Python backend (NumPy + SciPy) not available.",
          "Point RETICULATE_PYTHON at an environment that has them,",
          "or run diamondback_install_python().")
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

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
