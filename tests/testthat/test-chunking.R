# The reduction passes (stats, edges, overlaps) run in row chunks to bound
# memory. A patch spanning a chunk seam must measure exactly the same as one
# that does not. Real chunks are 4 million cells, so these tests shrink the
# chunk size to force seams through small rasters -- otherwise this code path
# would only ever be exercised on rasters too big to test.

with_chunk_cells <- function(n, code) {
  py <- diamondback:::db_py()
  old <- reticulate::py_to_r(py$CHUNK_CELLS)
  py$CHUNK_CELLS <- as.integer(n)
  on.exit(py$CHUNK_CELLS <- as.integer(old), add = TRUE)
  force(code)
}

test_that("metrics are identical whether or not patches cross chunk seams", {
  skip_if_no_python()
  set.seed(81)
  m <- matrix(rbinom(600, 1, 0.45), 30, 20)
  r <- rast_from(m)

  big <- patch_metrics(label_patches(r, quiet = TRUE), quiet = TRUE)$metrics
  # 40 cells per chunk over a 20-column raster = 2 rows per chunk, so nearly
  # every patch straddles a seam.
  small <- with_chunk_cells(40, {
    patch_metrics(label_patches(r, quiet = TRUE), quiet = TRUE)$metrics
  })

  expect_equal(big, small)
})

test_that("a single patch spanning many chunks is measured as one patch", {
  skip_if_no_python()
  # A vertical bar down the middle of a tall raster: one patch, crossing every
  # seam when chunks are two rows tall.
  m <- matrix(0, 30, 5); m[, 3] <- 1
  r <- rast_from(m)

  res <- with_chunk_cells(10, patch_metrics(label_patches(r, quiet = TRUE), quiet = TRUE))
  mt <- res$metrics

  expect_equal(nrow(mt), 1L)
  expect_equal(mt$cells, 30)
  # 30 cells in a 1-wide column: 30 edges on each side (60 vertical) plus the
  # top and bottom (2 horizontal) = 62 cell edges x 10 m.
  expect_equal(mt$edge_ns_cells, 60)
  expect_equal(mt$edge_ew_cells, 2)
  expect_equal(mt$perimeter, 620)
})

test_that("centroids are correct across chunk seams", {
  skip_if_no_python()
  m <- matrix(0, 30, 5); m[, 3] <- 1
  r <- rast_from(m)
  mt <- with_chunk_cells(10,
    patch_metrics(label_patches(r, quiet = TRUE), quiet = TRUE))$metrics
  # Column 3 of 5, rows 1-30 of a 30x5 grid with 10 m cells.
  expect_equal(mt$x_centroid, 25)
  expect_equal(mt$y_centroid, 150)
})

test_that("overlap cross-tabulation is identical across chunk sizes", {
  skip_if_no_python()
  set.seed(82)
  m1 <- rast_from(matrix(rbinom(600, 1, 0.45), 30, 20))
  m2 <- rast_from(matrix(rbinom(600, 1, 0.45), 30, 20))

  big <- compare_patches(m1, m2, quiet = TRUE)$overlaps
  small <- with_chunk_cells(40, compare_patches(m1, m2, quiet = TRUE)$overlaps)

  expect_equal(big, small)
  expect_gt(nrow(big), 5)   # the test is only meaningful if there are overlaps
})

test_that("lon/lat per-row weighting survives chunking", {
  skip_if_no_python()
  g <- terra::rast(nrows = 30, ncols = 20, xmin = 0, xmax = 20,
                   ymin = 20, ymax = 50, crs = "EPSG:4326")
  set.seed(83)
  terra::values(g) <- rbinom(600, 1, 0.45)

  big <- patch_metrics(label_patches(g, quiet = TRUE), quiet = TRUE)$metrics
  small <- with_chunk_cells(40,
    patch_metrics(label_patches(g, quiet = TRUE), quiet = TRUE))$metrics

  # Per-row edge lengths and cell areas are indexed by absolute row, so an
  # off-by-one at a chunk boundary would show up here as a latitude shift.
  expect_equal(big$area_m2, small$area_m2)
  expect_equal(big$perimeter, small$perimeter)
})

test_that("blockwise raster ingestion is independent of block size", {
  skip_if_no_python()
  set.seed(84)
  m <- matrix(rbinom(600, 1, 0.4), 30, 20)
  m[sample(600, 60)] <- NA
  r <- rast_from(m)

  ref <- label_patches(r, quiet = TRUE)

  # db_row_blocks drives reading from terra; force many small blocks.
  small_blocks <- testthat::with_mocked_bindings(
    label_patches(r, quiet = TRUE),
    db_row_blocks = function(nrow, ncol, target_cells = 4e6) {
      lapply(seq_len(nrow), function(s) list(row = s, nrows = 1L))
    },
    .package = "diamondback"
  )

  expect_equal(small_blocks$metadata$n_patches, ref$metadata$n_patches)
  expect_equal(small_blocks$metrics$cells, ref$metrics$cells)
  expect_equal(small_blocks$metadata$cells, ref$metadata$cells)
  expect_equal(terra::values(small_blocks$patches, mat = FALSE),
               terra::values(ref$patches, mat = FALSE))
})
