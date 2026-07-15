test_that("a single isolated patch is found", {
  skip_if_no_python()
  m <- matrix(0, 5, 5)
  m[2:3, 2:3] <- 1
  res <- label_patches(m, quiet = TRUE)
  expect_equal(res$metadata$n_patches, 1L)
  expect_equal(res$metrics$cells, 4)
})

test_that("multiple separate patches are found and numbered consecutively", {
  skip_if_no_python()
  m <- matrix(0, 5, 5)
  m[1, 1] <- 1; m[3, 3] <- 1; m[5, 5] <- 1
  res <- label_patches(m, quiet = TRUE)
  expect_equal(res$metadata$n_patches, 3L)
  expect_equal(sort(res$metrics$patch_id), 1:3)
  expect_equal(res$metrics$cells, c(1, 1, 1))
})

test_that("diagonal contact joins patches under 8 but not 4 connectivity", {
  skip_if_no_python()
  m <- matrix(c(1, 0,
                0, 1), nrow = 2, byrow = TRUE)
  expect_equal(label_patches(m, directions = 8, quiet = TRUE)$metadata$n_patches, 1L)
  expect_equal(label_patches(m, directions = 4, quiet = TRUE)$metadata$n_patches, 2L)
})

test_that("a diagonal-only chain is one patch under 8 and many under 4", {
  skip_if_no_python()
  m <- diag(5)
  expect_equal(label_patches(m, directions = 8, quiet = TRUE)$metadata$n_patches, 1L)
  expect_equal(label_patches(m, directions = 4, quiet = TRUE)$metadata$n_patches, 5L)
})

test_that("a hole inside a patch is background, not part of the patch", {
  skip_if_no_python()
  m <- matrix(1, 5, 5)
  m[3, 3] <- 0
  res <- label_patches(m, quiet = TRUE)
  expect_equal(res$metadata$n_patches, 1L)
  expect_equal(res$metrics$cells, 24)
  expect_equal(res$metadata$cells$background, 1)
})

test_that("an all-background raster gives zero patches and warns", {
  skip_if_no_python()
  m <- matrix(0, 4, 4)
  expect_warning(res <- label_patches(m, quiet = TRUE), "No foreground")
  expect_equal(res$metadata$n_patches, 0L)
  expect_equal(nrow(res$metrics), 0L)
})

test_that("an all-foreground raster is one patch covering every cell", {
  skip_if_no_python()
  m <- matrix(1, 4, 4)
  res <- label_patches(m, quiet = TRUE)
  expect_equal(res$metadata$n_patches, 1L)
  expect_equal(res$metrics$cells, 16)
  expect_equal(res$metadata$cells$background, 0)
})

test_that("a single-cell raster works in both states", {
  skip_if_no_python()
  expect_equal(label_patches(matrix(1, 1, 1), quiet = TRUE)$metadata$n_patches, 1L)
  expect_warning(
    z <- label_patches(matrix(0, 1, 1), quiet = TRUE)$metadata$n_patches
  )
  expect_equal(z, 0L)
})

# --- NA handling: the core promise of the package ---------------------------

test_that("NA is never foreground and never background by default", {
  skip_if_no_python()
  m <- matrix(c(1, 1, NA,
                0, 1, 0,
                NA, 0, 1), nrow = 3, byrow = TRUE)
  res <- label_patches(m, quiet = TRUE)
  cl <- res$metadata$cells
  expect_equal(cl$foreground, 4)
  expect_equal(cl$background, 3)
  expect_equal(cl$missing, 2)
  expect_equal(cl$outside, 0)
  expect_equal(cl$foreground + cl$background + cl$missing + cl$outside, 9)

  v <- terra::values(res$patches, mat = FALSE)
  expect_true(is.na(v[3]))   # NA stays NA
  expect_true(is.na(v[7]))
})

test_that("na = 'background' makes NA valid background, still never foreground", {
  skip_if_no_python()
  m <- matrix(c(1, NA, 0), nrow = 1)
  res <- label_patches(m, na = "background", quiet = TRUE)
  v <- terra::values(res$patches, mat = FALSE)
  expect_equal(v[2], 0)                       # NA became background
  expect_equal(res$metadata$n_patches, 1L)
  # The count of genuinely missing cells is still reported, not erased.
  expect_equal(res$metadata$cells$missing, 1)
})

test_that("an all-NA raster gives zero patches", {
  skip_if_no_python()
  m <- matrix(NA_real_, 3, 3)
  expect_warning(res <- label_patches(m, quiet = TRUE))
  expect_equal(res$metadata$n_patches, 0L)
  expect_equal(res$metadata$cells$missing, 9)
})

# --- masks ------------------------------------------------------------------

test_that("a mask excludes cells from the domain without making them background", {
  skip_if_no_python()
  v <- matrix(1, 3, 3)
  mk <- matrix(c(1, 1, 1,
                 1, 1, 1,
                 0, 0, 0), nrow = 3, byrow = TRUE)
  res <- label_patches(v, mask = mk, crop = FALSE, quiet = TRUE)
  expect_equal(res$metadata$cells$foreground, 6)
  expect_equal(res$metadata$cells$outside, 3)
  expect_equal(res$metadata$cells$background, 0)
  expect_equal(res$metadata$n_patches, 1L)
  expect_true(all(is.na(terra::values(res$patches, mat = FALSE)[7:9])))
})

test_that("an irregular mask splits what would otherwise be one patch", {
  skip_if_no_python()
  v <- matrix(1, 3, 3)
  mk <- matrix(c(1, 0, 1,
                 1, 0, 1,
                 1, 0, 1), nrow = 3, byrow = TRUE)
  res <- label_patches(v, mask = mk, crop = FALSE, directions = 8, quiet = TRUE)
  expect_equal(res$metadata$n_patches, 2L)
  expect_equal(res$metadata$cells$outside, 3)
})

test_that("NA in the mask means outside the domain", {
  skip_if_no_python()
  v <- matrix(1, 2, 2)
  mk <- matrix(c(1, NA, 1, 1), nrow = 2)
  res <- label_patches(v, mask = mk, crop = FALSE, quiet = TRUE)
  expect_equal(res$metadata$cells$outside, 1)
  expect_equal(res$metadata$cells$foreground, 3)
})

test_that("outside-domain beats missing when a cell is both", {
  skip_if_no_python()
  v <- matrix(c(NA, 1, 1, 1), nrow = 2)
  mk <- matrix(c(0, 1, 1, 1), nrow = 2)
  res <- label_patches(v, mask = mk, crop = FALSE, quiet = TRUE)
  expect_equal(res$metadata$cells$outside, 1)
  expect_equal(res$metadata$cells$missing, 0)
})

test_that("a misaligned mask errors rather than resampling silently", {
  skip_if_no_python()
  r <- rast_from(matrix(1, 4, 4))
  mk <- terra::rast(nrows = 4, ncols = 4, xmin = 5, xmax = 45, ymin = 0, ymax = 40,
                    crs = "EPSG:5070")
  terra::values(mk) <- 1
  expect_error(label_patches(r, mask = mk, quiet = TRUE), "does not align")
})

# --- classes ----------------------------------------------------------------

test_that("class selects one value from a categorical raster", {
  skip_if_no_python()
  m <- matrix(c(1, 1, 2,
                2, 1, 2,
                3, 3, 2), nrow = 3, byrow = TRUE)
  r1 <- label_patches(m, class = 1, quiet = TRUE)
  expect_equal(r1$metadata$cells$foreground, 3)
  expect_equal(r1$metadata$cells$background, 6)

  r2 <- label_patches(m, class = 2, quiet = TRUE)
  expect_equal(r2$metadata$cells$foreground, 4)
})

test_that("class = 'all' labels every value separately with unique IDs", {
  skip_if_no_python()
  m <- matrix(c(1, 1, 2,
                2, 1, 2,
                3, 3, 2), nrow = 3, byrow = TRUE)
  res <- label_patches(m, class = "all", directions = 4, quiet = TRUE)
  expect_true(all(c("patch_id", "class", "cells") %in% names(res$metrics)))
  expect_equal(sort(unique(res$metrics$class)), c("1", "2", "3"))
  expect_equal(sort(res$metrics$patch_id), seq_len(nrow(res$metrics)))
  expect_equal(sum(res$metrics$cells), 9)   # every cell belongs to some class
  expect_equal(res$metadata$cells$background, 0)
})

test_that("a class vector selects only the values requested, as one class", {
  skip_if_no_python()
  m <- matrix(c(1, 2, 3, 1), nrow = 2)
  res <- label_patches(m, class = c(1, 3), quiet = TRUE)
  expect_equal(sum(res$metrics$cells), 3)
  expect_equal(res$metadata$cells$background, 1)  # the 2 is background
  # c(1, 3) is one class, so there is no class column and the values it names
  # are labelled together where they touch.
  expect_false("class" %in% names(res$metrics))
})

test_that("NA as a class is rejected", {
  skip_if_no_python()
  expect_error(label_patches(matrix(1, 2, 2), class = NA, quiet = TRUE), "NA")
  expect_error(label_patches(matrix(1, 2, 2), class = c(1, NA), quiet = TRUE), "NA")
})

test_that("invalid directions are rejected", {
  skip_if_no_python()
  expect_error(label_patches(matrix(1, 2, 2), directions = 6, quiet = TRUE),
               "must be")
})

# --- outputs ----------------------------------------------------------------

test_that("output writes a raster and the result reads back from it", {
  skip_if_no_python()
  f <- withr::local_tempfile(fileext = ".tif")
  m <- matrix(c(1, 0, 0, 1), nrow = 2)
  res <- label_patches(m, output = f, quiet = TRUE)
  expect_true(file.exists(f))
  expect_equal(res$path, f)
  expect_s4_class(res$patches, "SpatRaster")
  expect_equal(terra::ncell(res$patches), 4)
})

test_that("output refuses to clobber an existing file unless told to", {
  skip_if_no_python()
  f <- withr::local_tempfile(fileext = ".tif")
  m <- matrix(c(1, 0, 0, 1), nrow = 2)
  label_patches(m, output = f, quiet = TRUE)
  expect_error(label_patches(m, output = f, quiet = TRUE), "already exists")
  expect_no_error(label_patches(m, output = f, overwrite = TRUE, quiet = TRUE))
})

test_that("cropping to a mask reduces the grid but not the patch count", {
  skip_if_no_python()
  v <- matrix(0, 20, 20); v[9:11, 9:11] <- 1
  mk <- matrix(0, 20, 20); mk[8:12, 8:12] <- 1
  full <- label_patches(v, mask = mk, crop = FALSE, quiet = TRUE)
  crop <- label_patches(v, mask = mk, crop = TRUE, quiet = TRUE)
  expect_equal(full$metadata$n_patches, crop$metadata$n_patches)
  expect_equal(full$metrics$cells, crop$metrics$cells)
  expect_lt(crop$metadata$geometry$ncell, full$metadata$geometry$ncell)
})
