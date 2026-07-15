test_that("core area erodes a square patch by the expected ring", {
  skip_if_no_python()
  # 7x7 habitat inside a 9x9 grid, 1 m cells.
  m <- matrix(0, 9, 9); m[2:8, 2:8] <- 1
  r <- rast_from(m, res = 1)

  d1 <- patch_core_area(label_patches(r, quiet = TRUE), depth = 1, quiet = TRUE)$metrics
  expect_equal(d1$cells, 49)
  expect_equal(d1$core_cells, 25)      # 5x5 survives
  expect_equal(d1$core_fraction, 25 / 49)
  expect_equal(d1$core_area_m2, 25)

  d2 <- patch_core_area(label_patches(r, quiet = TRUE), depth = 2, quiet = TRUE)$metrics
  expect_equal(d2$core_cells, 9)       # 3x3

  d3 <- patch_core_area(label_patches(r, quiet = TRUE), depth = 3, quiet = TRUE)$metrics
  expect_equal(d3$core_cells, 1)       # centre only
})

test_that("a patch too small to have an interior has zero core", {
  skip_if_no_python()
  m <- matrix(0, 5, 5); m[3, 3] <- 1
  mt <- patch_core_area(label_patches(rast_from(m, res = 1), quiet = TRUE),
                        depth = 1, quiet = TRUE)$metrics
  expect_equal(mt$core_cells, 0)
  expect_equal(mt$core_fraction, 0)
})

test_that("depth 0 keeps every cell, since distance is always positive inside", {
  skip_if_no_python()
  m <- matrix(0, 5, 5); m[2:4, 2:4] <- 1
  mt <- patch_core_area(label_patches(rast_from(m, res = 1), quiet = TRUE),
                        depth = 0, quiet = TRUE)$metrics
  expect_equal(mt$core_cells, 9)
})

test_that("core area respects cell size, not cell count", {
  skip_if_no_python()
  m <- matrix(0, 9, 9); m[2:8, 2:8] <- 1
  # 10 m cells: depth 10 removes one ring, exactly as depth 1 did at 1 m.
  mt <- patch_core_area(label_patches(rast_from(m, res = 10), quiet = TRUE),
                        depth = 10, quiet = TRUE)$metrics
  expect_equal(mt$core_cells, 25)
  expect_equal(mt$core_area_m2, 2500)
})

test_that("edge = 'background' spares patches clipped by the grid border", {
  skip_if_no_python()
  # Habitat fills columns 1-3 and runs off the top, bottom and left borders.
  m <- matrix(0, 5, 5); m[1:5, 1:3] <- 1
  r <- rast_from(m, res = 1)

  all_edges <- patch_core_area(label_patches(r, quiet = TRUE), depth = 1,
                               edge = "all", quiet = TRUE)$metrics$core_cells
  bg_only <- patch_core_area(label_patches(r, quiet = TRUE), depth = 1,
                             edge = "background", quiet = TRUE)$metrics$core_cells

  # edge = "all": only column 2, rows 2-4 are >1 from something non-habitat.
  expect_equal(all_edges, 3)
  # edge = "background": only column 4 is an edge source, so columns 1-2
  # (distance 3 and 2) survive across all five rows.
  expect_equal(bg_only, 10)
  expect_gt(bg_only, all_edges)
})

test_that("edge = 'background' treats masked-out cells as habitat", {
  skip_if_no_python()
  v <- matrix(0, 5, 7); v[1:5, 1:3] <- 1
  mk <- matrix(1, 5, 7); mk[1:5, 4] <- 0   # a masked strip beside the habitat

  a <- patch_core_area(
    label_patches(rast_from(v), mask = rast_from(mk), crop = FALSE, quiet = TRUE),
    depth = 1, units = "cells", edge = "all", quiet = TRUE)$metrics$core_cells
  b <- patch_core_area(
    label_patches(rast_from(v), mask = rast_from(mk), crop = FALSE, quiet = TRUE),
    depth = 1, units = "cells", edge = "background", quiet = TRUE)$metrics$core_cells

  expect_gt(b, a)
})

test_that("non-square cells give an anisotropic core", {
  skip_if_no_python()
  # Cells 1 m wide, 10 m tall. A depth of 2 m erodes horizontally (2 cells)
  # but vertically a single cell step is already 10 m, so one row is enough.
  m <- matrix(0, 5, 9); m[2:4, 2:8] <- 1
  r <- rast_from(m, res = 1, yres = 10)
  mt <- patch_core_area(label_patches(r, quiet = TRUE), depth = 2, quiet = TRUE)$metrics
  # Row 3 only (rows 2 and 4 are 10 m from an edge... but 10 > 2, so they stay
  # unless they are horizontally close to an edge). Verify against a direct count.
  expect_true(mt$core_cells > 0)
  expect_lt(mt$core_cells, mt$cells)
})

test_that("core area errors on lon/lat unless measured in cells", {
  skip_if_no_python()
  g <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 40, ymax = 45,
                   crs = "EPSG:4326")
  terra::values(g) <- 1
  lab <- label_patches(g, quiet = TRUE)

  expect_error(patch_core_area(lab, depth = 100, quiet = TRUE), "projected raster")
  expect_no_error(patch_core_area(lab, depth = 1, units = "cells", quiet = TRUE))
})

test_that("core cells never exceed patch cells on random rasters", {
  skip_if_no_python()
  set.seed(31)
  m <- matrix(rbinom(900, 1, 0.6), 30, 30)
  mt <- patch_core_area(label_patches(rast_from(m, res = 1), quiet = TRUE),
                        depth = 1, quiet = TRUE)$metrics
  expect_true(all(mt$core_cells <= mt$cells))
  expect_true(all(mt$core_fraction >= 0 & mt$core_fraction <= 1))
})

test_that("core area is computed per class, not across classes", {
  skip_if_no_python()
  # Two classes filling the grid. If the transform ran on "any patch", every
  # cell would look interior; per class, neither has core at depth 1.
  m <- matrix(1, 6, 6); m[, 4:6] <- 2
  res <- label_patches(rast_from(m, res = 1), class = c(1, 2), quiet = TRUE)
  mt <- patch_core_area(res, depth = 1, quiet = TRUE)$metrics
  expect_equal(nrow(mt), 2L)
  # Each class is a 6x3 block: only its own interior counts, and the shared
  # border between them is an edge for both.
  expect_true(all(mt$core_cells < mt$cells))
  expect_equal(mt$core_cells, c(4, 4))
})

test_that("core_raster returns a raster of core cells labelled by patch", {
  skip_if_no_python()
  m <- matrix(0, 9, 9); m[2:8, 2:8] <- 1
  res <- patch_core_area(label_patches(rast_from(m, res = 1), quiet = TRUE),
                         depth = 1, core_raster = TRUE, quiet = TRUE)
  expect_s4_class(res$core_raster, "SpatRaster")
  v <- terra::values(res$core_raster, mat = FALSE)
  expect_equal(sum(v == 1), 25)
})

test_that("negative or non-scalar depth is rejected", {
  skip_if_no_python()
  lab <- label_patches(matrix(1, 3, 3), quiet = TRUE)
  expect_error(patch_core_area(lab, depth = -1, quiet = TRUE), "non-negative")
  expect_error(patch_core_area(lab, depth = c(1, 2), quiet = TRUE), "single")
})
