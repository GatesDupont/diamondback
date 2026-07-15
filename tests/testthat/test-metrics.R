# Hand-calculated answers. Every expectation here is derivable with a pencil.

test_that("area, perimeter and edges are exact for an interior square", {
  skip_if_no_python()
  # 2x2 patch surrounded by background, 10 m cells.
  m <- matrix(0, 4, 4); m[2:3, 2:3] <- 1
  res <- patch_metrics(label_patches(rast_from(m), quiet = TRUE), quiet = TRUE)
  mt <- res$metrics

  expect_equal(mt$cells, 4)
  expect_equal(mt$area_m2, 400)          # 4 cells x 100 m2
  expect_equal(mt$area_ha, 0.04)
  expect_equal(mt$area_km2, 4e-04)
  expect_equal(mt$perimeter, 80)         # 8 exposed cell edges x 10 m
  expect_equal(mt$edge_valid, 80)        # all against real background
  expect_equal(mt$edge_domain, 0)
  expect_false(mt$touches_domain_edge)
  expect_equal(mt$edge_ns_cells, 4)
  expect_equal(mt$edge_ew_cells, 4)
  expect_equal(mt$edge_area_ratio, 80 / 400)
})

test_that("a single cell has a perimeter of four cell edges", {
  skip_if_no_python()
  m <- matrix(0, 3, 3); m[2, 2] <- 1
  mt <- patch_metrics(label_patches(rast_from(m), quiet = TRUE), quiet = TRUE)$metrics
  expect_equal(mt$cells, 1)
  expect_equal(mt$perimeter, 40)
  expect_equal(mt$edge_ns_cells, 2)
  expect_equal(mt$edge_ew_cells, 2)
})

test_that("the grid border counts as domain edge, not habitat edge", {
  skip_if_no_python()
  m <- matrix(0, 4, 4); m[1:2, 1:2] <- 1   # corner patch
  mt <- patch_metrics(label_patches(rast_from(m), quiet = TRUE), quiet = TRUE)$metrics
  expect_equal(mt$perimeter, 80)
  expect_equal(mt$edge_valid, 40)    # 2 sides face background
  expect_equal(mt$edge_domain, 40)   # 2 sides face off-grid
  expect_true(mt$touches_domain_edge)
})

test_that("edge against a masked-out cell is domain edge, not habitat edge", {
  skip_if_no_python()
  v <- matrix(c(1, 0, 0,
                0, 0, 0,
                0, 0, 0), nrow = 3, byrow = TRUE)
  mk <- matrix(1, 3, 3); mk[1, 2] <- 0   # mask out the cell to the patch's right
  res <- patch_metrics(label_patches(rast_from(v), mask = rast_from(mk),
                                     crop = FALSE, quiet = TRUE), quiet = TRUE)
  mt <- res$metrics
  # Patch is at the top-left corner: top and left are grid border (domain),
  # right is masked out (domain), bottom is background (valid).
  expect_equal(mt$perimeter, 40)
  expect_equal(mt$edge_valid, 10)
  expect_equal(mt$edge_domain, 30)
})

test_that("edge against a missing cell is domain edge by default and habitat edge with na='background'", {
  skip_if_no_python()
  v <- matrix(c(1, NA, 0,
                0, 0, 0,
                0, 0, 0), nrow = 3, byrow = TRUE)
  a <- patch_metrics(label_patches(rast_from(v), quiet = TRUE), quiet = TRUE)$metrics
  expect_equal(a$edge_valid, 10)     # only the cell below
  expect_equal(a$edge_domain, 30)    # top, left border + the NA to the right

  b <- patch_metrics(label_patches(rast_from(v), na = "background", quiet = TRUE),
                     quiet = TRUE)$metrics
  expect_equal(b$edge_valid, 20)     # the NA now counts as habitat edge
  expect_equal(b$edge_domain, 20)    # only the grid border remains
})

test_that("non-square cells use the right edge length on each axis", {
  skip_if_no_python()
  # Cells 10 m wide (E-W) and 2 m tall (N-S). One isolated cell:
  # 2 vertical edges of 2 m + 2 horizontal edges of 10 m = 24 m.
  m <- matrix(0, 3, 3); m[2, 2] <- 1
  r <- rast_from(m, res = 10, yres = 2)
  mt <- patch_metrics(label_patches(r, quiet = TRUE), quiet = TRUE)$metrics
  expect_equal(mt$area_m2, 20)
  expect_equal(mt$perimeter, 24)
  expect_equal(mt$edge_ns_cells, 2)
  expect_equal(mt$edge_ew_cells, 2)
})

test_that("non-square cells: a wide patch is dominated by its long edges", {
  skip_if_no_python()
  # 1x3 horizontal strip, cells 10 wide x 2 tall.
  # vertical edges: 2 (ends) x 2 m = 4; horizontal: 6 x 10 m = 60. Total 64.
  m <- matrix(0, 3, 5); m[2, 2:4] <- 1
  r <- rast_from(m, res = 10, yres = 2)
  mt <- patch_metrics(label_patches(r, quiet = TRUE), quiet = TRUE)$metrics
  expect_equal(mt$perimeter, 64)
  expect_equal(mt$edge_ns_cells, 2)
  expect_equal(mt$edge_ew_cells, 6)
})

test_that("a patch with a hole has the hole's edge in its perimeter", {
  skip_if_no_python()
  m <- matrix(1, 5, 5); m[3, 3] <- 0
  mt <- patch_metrics(label_patches(rast_from(m), quiet = TRUE), quiet = TRUE)$metrics
  expect_equal(mt$cells, 24)
  # outer boundary: 20 cell edges (all grid border) + hole: 4 = 24 edges
  expect_equal(mt$perimeter, 240)
  expect_equal(mt$edge_valid, 40)     # only the hole is real habitat edge
  expect_equal(mt$edge_domain, 200)   # the outer ring is all grid border
})

test_that("bounding box and centroid are in map coordinates", {
  skip_if_no_python()
  m <- matrix(0, 4, 4); m[2:3, 2:3] <- 1
  mt <- patch_metrics(label_patches(rast_from(m), quiet = TRUE), quiet = TRUE)$metrics
  # Rows 2-3, cols 2-3 of a 4x4 grid with 10 m cells, ymax = 40.
  expect_equal(mt$xmin, 10); expect_equal(mt$xmax, 30)
  expect_equal(mt$ymin, 10); expect_equal(mt$ymax, 30)
  expect_equal(mt$x_centroid, 20)
  expect_equal(mt$y_centroid, 20)
})

test_that("the centroid of an L-shaped patch may lie outside it", {
  skip_if_no_python()
  m <- matrix(0, 4, 4)
  m[2, 2] <- 1; m[3, 2] <- 1; m[3, 3] <- 1
  mt <- patch_metrics(label_patches(rast_from(m), quiet = TRUE), quiet = TRUE)$metrics
  # mean col index = (1+1+2)/3, mean row = (1+2+2)/3 -> centroid at (~23.3, ~23.3)
  expect_equal(mt$x_centroid, 10 * (4 / 3) + 5, tolerance = 1e-8)
  expect_equal(mt$y_centroid, 40 - 10 * (5 / 3) - 5, tolerance = 1e-8)
})

test_that("units = 'cells' reports counts rather than map units", {
  skip_if_no_python()
  m <- matrix(0, 4, 4); m[2:3, 2:3] <- 1
  mt <- patch_metrics(label_patches(rast_from(m), quiet = TRUE),
                      units = "cells", quiet = TRUE)$metrics
  expect_equal(mt$area_cells, 4)
  expect_equal(mt$perimeter, 8)
  expect_false("area_m2" %in% names(mt))
})

test_that("units = 'km' scales lengths and areas consistently", {
  skip_if_no_python()
  m <- matrix(0, 4, 4); m[2:3, 2:3] <- 1
  mt <- patch_metrics(label_patches(rast_from(m), quiet = TRUE),
                      units = "km", quiet = TRUE)$metrics
  expect_equal(mt$perimeter, 0.08)      # 80 m
  expect_equal(mt$area_km2, 4e-04)
})

test_that("edge components always sum to the perimeter", {
  skip_if_no_python()
  set.seed(21)
  m <- matrix(rbinom(400, 1, 0.4), 20, 20)
  mt <- patch_metrics(label_patches(rast_from(m), quiet = TRUE), quiet = TRUE)$metrics
  expect_equal(mt$edge_valid + mt$edge_domain, mt$perimeter)
  # And cell counts are consistent with lengths on a square grid.
  expect_equal((mt$edge_ns_cells + mt$edge_ew_cells) * 10, mt$perimeter)
})

test_that("summed patch cells equal the foreground count", {
  skip_if_no_python()
  set.seed(22)
  m <- matrix(rbinom(400, 1, 0.35), 20, 20)
  res <- patch_metrics(label_patches(rast_from(m), quiet = TRUE), quiet = TRUE)
  expect_equal(sum(res$metrics$cells), res$metadata$cells$foreground)
})

test_that("an empty result gives a zero-row metrics table with the right columns", {
  skip_if_no_python()
  expect_warning(lab <- label_patches(matrix(0, 4, 4), quiet = TRUE))
  mt <- patch_metrics(lab, quiet = TRUE)$metrics
  expect_equal(nrow(mt), 0L)
  expect_true(all(c("patch_id", "cells", "perimeter", "touches_domain_edge") %in% names(mt)))
})

# --- geographic rasters -----------------------------------------------------

test_that("lon/lat area matches terra::cellSize exactly", {
  skip_if_no_python()
  g <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10,
                   ymin = 30, ymax = 40, crs = "EPSG:4326")
  set.seed(1)
  terra::values(g) <- rbinom(100, 1, 0.5)

  res <- patch_metrics(label_patches(g, quiet = TRUE), quiet = TRUE)
  cs <- terra::cellSize(g, unit = "m")
  expected <- sum(terra::values(cs, mat = FALSE)[terra::values(g, mat = FALSE) == 1])
  expect_equal(sum(res$metrics$area_m2), expected, tolerance = 1e-6)
})

test_that("lon/lat perimeter uses geodesic per-row edge lengths", {
  skip_if_no_python()
  g <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5,
                   ymin = 40, ymax = 45, crs = "EPSG:4326")
  terra::values(g) <- 0
  g[3, 3] <- 1   # one cell spanning lon 2-3, lat 42-43

  mt <- patch_metrics(label_patches(g, quiet = TRUE), quiet = TRUE)$metrics

  # Independent expectation: two N-S edges plus the E-W widths at both the
  # top (43) and bottom (42) latitudes of the cell.
  dy <- terra::distance(cbind(0, 43), cbind(0, 42), lonlat = TRUE, pairwise = TRUE)
  dx_top <- terra::distance(cbind(0, 43), cbind(1, 43), lonlat = TRUE, pairwise = TRUE)
  dx_bot <- terra::distance(cbind(0, 42), cbind(1, 42), lonlat = TRUE, pairwise = TRUE)
  expect_equal(mt$perimeter, as.numeric(2 * dy + dx_top + dx_bot), tolerance = 1e-6)
})

test_that("lon/lat cell area varies with latitude, and diamondback follows it", {
  skip_if_no_python()
  g <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 30,
                   ymin = 0, ymax = 30, crs = "EPSG:4326")
  terra::values(g) <- 0
  g[1, 2] <- 1   # high-latitude cell
  north <- patch_metrics(label_patches(g, quiet = TRUE), quiet = TRUE)$metrics$area_m2

  terra::values(g) <- 0
  g[2, 2] <- 1   # equatorial cell
  south <- patch_metrics(label_patches(g, quiet = TRUE), quiet = TRUE)$metrics$area_m2

  expect_lt(north, south)   # cells shrink towards the pole
})
