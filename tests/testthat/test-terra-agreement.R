# The whole premise of the package is that SciPy gives the same answer as
# terra::patches(), only faster. If that is not true, nothing else matters.

test_that("patch counts match terra::patches() on random rasters", {
  skip_if_no_python()
  for (dir in c(4, 8)) {
    for (seed in 1:5) {
      set.seed(seed)
      m <- matrix(rbinom(400, 1, 0.45), 20, 20)
      db <- label_patches(m, directions = dir, quiet = TRUE)$metadata$n_patches
      tr <- terra_patch_count(m, directions = dir)
      expect_equal(db, tr,
                   info = sprintf("directions = %d, seed = %d", dir, seed))
    }
  }
})

test_that("the partition of cells matches terra::patches(), not just the count", {
  skip_if_no_python()
  set.seed(99)
  m <- matrix(rbinom(400, 1, 0.4), 20, 20)

  for (dir in c(4, 8)) {
    db <- label_patches(m, directions = dir, quiet = TRUE)
    db_v <- terra::values(db$patches, mat = FALSE)
    db_v[is.na(db_v)] <- 0

    r <- rast_from(m)
    r[r == 0] <- NA
    tr_v <- terra::values(terra::patches(r, directions = dir, zeroAsNA = FALSE),
                          mat = FALSE)

    expect_true(same_partition(db_v, tr_v),
                info = sprintf("directions = %d", dir))
  }
})

test_that("agreement holds at different foreground densities", {
  skip_if_no_python()
  for (p in c(0.05, 0.2, 0.5, 0.8, 0.95)) {
    set.seed(round(p * 100))
    m <- matrix(rbinom(900, 1, p), 30, 30)
    expect_equal(
      label_patches(m, directions = 8, quiet = TRUE)$metadata$n_patches,
      terra_patch_count(m, directions = 8),
      info = sprintf("density = %.2f", p)
    )
  }
})

test_that("agreement holds when NA cells are present", {
  skip_if_no_python()
  set.seed(11)
  m <- matrix(rbinom(400, 1, 0.4), 20, 20)
  m[sample(400, 40)] <- NA

  db <- label_patches(m, directions = 8, quiet = TRUE)$metadata$n_patches

  # terra: NA is already no-data, and 0 must become no-data too, which leaves
  # exactly diamondback's default reading of the same raster.
  r <- rast_from(m)
  r[r == 0] <- NA
  tr <- length(unique(stats::na.omit(
    terra::values(terra::patches(r, directions = 8, zeroAsNA = FALSE), mat = FALSE))))

  expect_equal(db, tr)
})

test_that("patch cell counts match terra::patches() patch by patch", {
  skip_if_no_python()
  set.seed(5)
  m <- matrix(rbinom(400, 1, 0.4), 20, 20)

  db <- label_patches(m, directions = 8, quiet = TRUE)
  r <- rast_from(m)
  r[r == 0] <- NA
  tr <- terra::patches(r, directions = 8, zeroAsNA = FALSE)

  # Labels differ in numbering, so compare the sorted size distributions.
  db_sizes <- sort(db$metrics$cells)
  tr_sizes <- sort(as.numeric(table(terra::values(tr, mat = FALSE))))
  expect_equal(db_sizes, tr_sizes)
})
