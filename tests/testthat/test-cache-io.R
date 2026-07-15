make_src <- function(seed = 1, n = 20) {
  f <- tempfile(fileext = ".tif")
  r <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n * 10,
                   ymin = 0, ymax = n * 10, crs = "EPSG:5070")
  set.seed(seed)
  terra::values(r) <- rbinom(n * n, 1, 0.4)
  terra::writeRaster(r, f, overwrite = TRUE)
  f
}

test_that("write_patch_result and read_patch_result round-trip", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  m <- matrix(c(1, 1, 0, 0, 1, 0, 0, 1, 1), nrow = 3)
  res <- analyze_patches(m, quiet = TRUE)

  write_patch_result(res, d, overwrite = TRUE, quiet = TRUE)
  expect_true(file.exists(file.path(d, "patches.tif")))
  expect_true(file.exists(file.path(d, "metrics.csv")))
  expect_true(file.exists(file.path(d, "metadata.json")))

  back <- read_patch_result(d)
  expect_s3_class(back, "patch_result")
  expect_equal(back$metadata$n_patches, res$metadata$n_patches)
  expect_equal(back$metrics$cells, res$metrics$cells)
  expect_equal(terra::values(back$patches, mat = FALSE),
               terra::values(res$patches, mat = FALSE))
})

test_that("write_patch_result refuses to clobber without overwrite", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  res <- analyze_patches(matrix(c(1, 0, 0, 1), 2, 2), quiet = TRUE)
  write_patch_result(res, d, overwrite = TRUE, quiet = TRUE)
  expect_error(write_patch_result(res, d, quiet = TRUE), "already contains")
})

test_that("read_patch_result errors clearly on an incomplete directory", {
  d <- withr::local_tempdir()
  expect_error(read_patch_result(d), "not a complete patch result")
})

test_that("a restored result still supports metrics and comparison", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  m <- matrix(0, 6, 6); m[2:4, 2:4] <- 1
  write_patch_result(analyze_patches(rast_from(m), quiet = TRUE), d,
                     overwrite = TRUE, quiet = TRUE)

  # The NumPy arrays are gone; they must be rebuilt from the raster.
  back <- read_patch_result(d)
  again <- patch_metrics(back, quiet = TRUE)
  expect_equal(again$metrics$cells, 9)
  expect_equal(again$metrics$perimeter, 120)

  cmp <- compare_patches(back, back, quiet = TRUE)
  expect_equal(cmp$overlaps$overlap_cells, 9)
})

# --- caching ----------------------------------------------------------------

test_that("an identical request hits the cache", {
  skip_if_no_python()
  f <- make_src(); d <- withr::local_tempdir()
  first <- analyze_patches(f, output_dir = d, quiet = TRUE)
  expect_message(second <- analyze_patches(f, output_dir = d, quiet = FALSE),
                 "Reusing cached result")
  expect_equal(second$metadata$n_patches, first$metadata$n_patches)
  expect_equal(second$metrics$cells, first$metrics$cells)
})

test_that("changing connectivity misses the cache and says so", {
  skip_if_no_python()
  f <- make_src(); d <- withr::local_tempdir()
  analyze_patches(f, directions = 8, output_dir = d, quiet = TRUE)
  expect_message(analyze_patches(f, directions = 4, output_dir = d, quiet = FALSE),
                 "directions differs")
})

test_that("changing the class misses the cache", {
  skip_if_no_python()
  f <- make_src(); d <- withr::local_tempdir()
  analyze_patches(f, output_dir = d, quiet = TRUE)
  expect_message(analyze_patches(f, class = 1, output_dir = d, quiet = FALSE),
                 "class differs")
})

test_that("changing the NA policy misses the cache", {
  skip_if_no_python()
  f <- make_src(); d <- withr::local_tempdir()
  analyze_patches(f, output_dir = d, quiet = TRUE)
  expect_message(analyze_patches(f, na = "background", output_dir = d, quiet = FALSE),
                 "na differs")
})

test_that("a modified source file misses the cache", {
  skip_if_no_python()
  f <- make_src(seed = 1); d <- withr::local_tempdir()
  analyze_patches(f, output_dir = d, quiet = TRUE)

  r <- terra::rast(f)
  set.seed(2); terra::values(r) <- rbinom(terra::ncell(r), 1, 0.7)
  terra::writeRaster(r, f, overwrite = TRUE)

  expect_message(analyze_patches(f, output_dir = d, quiet = FALSE), "source")
})

test_that("adding a mask misses the cache", {
  skip_if_no_python()
  f <- make_src(); d <- withr::local_tempdir()
  analyze_patches(f, output_dir = d, quiet = TRUE)

  mk <- terra::rast(f); terra::values(mk) <- 1
  fm <- tempfile(fileext = ".tif"); terra::writeRaster(mk, fm, overwrite = TRUE)
  expect_message(analyze_patches(f, mask = fm, output_dir = d, quiet = FALSE),
                 "mask_source")
})

test_that("cache = FALSE always recomputes", {
  skip_if_no_python()
  f <- make_src(); d <- withr::local_tempdir()
  analyze_patches(f, output_dir = d, quiet = TRUE)
  expect_no_message(analyze_patches(f, output_dir = d, cache = FALSE, quiet = TRUE))
})

test_that("an in-memory raster is never served from cache", {
  skip_if_no_python()
  # No file behind it means no way to prove it did not change, so a hit would
  # be a guess. The package declines to guess.
  d <- withr::local_tempdir()
  r <- terra::rast(nrows = 1200, ncols = 1200, xmin = 0, xmax = 1200,
                   ymin = 0, ymax = 1200, crs = "EPSG:5070")
  set.seed(4); terra::values(r) <- rbinom(terra::ncell(r), 1, 0.3)

  analyze_patches(r, output_dir = d, metrics = FALSE, quiet = TRUE)
  expect_message(analyze_patches(r, output_dir = d, metrics = FALSE, quiet = FALSE),
                 "cannot be fingerprinted")
})

test_that("a small in-memory raster is fingerprinted by content and can hit", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  r <- rast_from(matrix(c(1, 0, 0, 1), 2, 2))
  analyze_patches(r, output_dir = d, quiet = TRUE)
  expect_message(analyze_patches(r, output_dir = d, quiet = FALSE),
                 "Reusing cached result")
})

test_that("the cache key ignores things that cannot change the labels", {
  skip_if_no_python()
  f <- make_src(); d <- withr::local_tempdir()
  analyze_patches(f, output_dir = d, units = "m", quiet = TRUE)
  # Different metric units, same labels: this must still hit.
  expect_message(analyze_patches(f, output_dir = d, units = "km", quiet = FALSE),
                 "Reusing cached result")
})

test_that("metadata records what is needed to reproduce a run", {
  skip_if_no_python()
  f <- make_src(); d <- withr::local_tempdir()
  res <- analyze_patches(f, class = 1, directions = 4, output_dir = d, quiet = TRUE)
  m <- read_patch_result(d)$metadata

  expect_equal(m$directions, 4)
  # class is a group specification now, not a bare value.
  expect_equal(m$class$labels, "1")
  expect_equal(unname(unlist(m$class$groups)), 1)
  expect_equal(m$na, "outside")
  expect_equal(m$algorithm_version, "2")
  expect_false(is.null(m$source$hash))
  expect_false(is.null(m$geometry$crs))
  expect_false(is.null(m$backend$scipy))
  expect_false(is.null(m$package_version))
  expect_equal(m$cells$foreground + m$cells$background + m$cells$missing +
                 m$cells$outside, m$geometry$ncell)
})
