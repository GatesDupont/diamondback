test_that("SpatRaster, filename and matrix inputs agree", {
  skip_if_no_python()
  m <- matrix(c(1, 1, 0, 0, 1, 0, 0, 1, 1), nrow = 3)
  r <- rast_from(m)
  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)

  a <- label_patches(m, quiet = TRUE)$metadata$n_patches
  b <- label_patches(r, quiet = TRUE)$metadata$n_patches
  d <- label_patches(f, quiet = TRUE)$metadata$n_patches
  expect_equal(a, b)
  expect_equal(b, d)
})

test_that("a missing file is reported as a missing file", {
  expect_error(label_patches("definitely_not_here.tif", quiet = TRUE),
               "does not exist")
})

test_that("an unsupported input type is rejected clearly", {
  expect_error(label_patches(list(1, 2), quiet = TRUE), "must be a")
  expect_error(label_patches(1:10, quiet = TRUE), "must be a")
})

test_that("a multi-layer raster is rejected with a useful suggestion", {
  skip_if_no_python()
  st <- c(rast_from(matrix(1, 2, 2)), rast_from(matrix(1, 2, 2)))
  expect_error(label_patches(st, quiet = TRUE), "one layer at a time")
})

test_that("passing a patch_result where a raster is expected says what to do", {
  skip_if_no_python()
  res <- analyze_patches(matrix(1, 2, 2), quiet = TRUE)
  expect_error(db_as_rast(res), "\\$patches")
})

# --- geometry comparison ----------------------------------------------------

test_that("origin shifts are detected even at identical resolution", {
  g1 <- db_geometry(terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 40,
                                ymin = 0, ymax = 40, crs = "EPSG:5070"))
  g2 <- db_geometry(terra::rast(nrows = 4, ncols = 4, xmin = 5, xmax = 45,
                                ymin = 0, ymax = 40, crs = "EPSG:5070"))
  expect_true("extent" %in% db_geometry_diff(g1, g2))
})

test_that("identical geometries compare equal", {
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 40, ymin = 0, ymax = 40,
                   crs = "EPSG:5070")
  expect_length(db_geometry_diff(db_geometry(r), db_geometry(r)), 0L)
})

test_that("different CRS is detected", {
  g1 <- db_geometry(terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 40,
                                ymin = 0, ymax = 40, crs = "EPSG:5070"))
  g2 <- db_geometry(terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 40,
                                ymin = 0, ymax = 40, crs = "EPSG:3857"))
  expect_true("CRS" %in% db_geometry_diff(g1, g2))
})

test_that("an unknown CRS is treated as planar rather than guessed from extent", {
  # A bare matrix has a 0-4 extent, which perhaps=TRUE would read as degrees.
  r <- db_as_rast(matrix(1, 4, 4))
  expect_false(db_is_lonlat(r))
})

# --- memory -----------------------------------------------------------------

test_that("memory estimates scale with cell count and stage", {
  n <- 1e6
  expect_equal(db_estimate_memory(n, "label"), n * 6)
  expect_gt(db_estimate_memory(n, "core"), db_estimate_memory(n, "label"))
})

test_that("an oversized allocation errors before any work happens", {
  expect_error(
    db_check_memory(1e12, "label", max_memory_frac = 0.6,
                    memory_limit = 1e9, quiet = TRUE),
    class = "diamondback_memory_error"
  )
})

test_that("the memory error names the estimate and the ways out", {
  err <- tryCatch(
    db_check_memory(1e12, "label", memory_limit = 1e9, quiet = TRUE),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "exceeds the safe limit")
  expect_match(err, "max_memory_frac")
  expect_match(err, "memory_limit")
})

test_that("label_patches respects an explicit memory limit", {
  skip_if_no_python()
  expect_error(
    label_patches(matrix(1, 100, 100), memory_limit = 1000, quiet = TRUE),
    class = "diamondback_memory_error"
  )
})

test_that("a comfortable allocation passes", {
  expect_silent(db_check_memory(100, "label", memory_limit = 1e9, quiet = TRUE))
})

test_that("db_memory_report returns per-stage estimates", {
  r <- rast_from(matrix(1, 10, 10))
  est <- suppressMessages(db_memory_report(r))
  expect_named(est, c("label", "metrics", "core"))
  expect_equal(est[["label"]], 100 * 6)
})

# --- environment ------------------------------------------------------------

test_that("diamondback_check reports a working environment", {
  skip_if_no_python()
  chk <- diamondback_check(verbose = FALSE)
  expect_true(chk$ok)
  expect_true(chk$test_passed)
  expect_false(is.na(chk$numpy))
  expect_false(is.na(chk$scipy))
  expect_match(chk$r, "^4\\.")
})

test_that("the temp-file registry only removes files diamondback made", {
  f_mine <- file.path(tempdir(), "diamondback_registry_test.tif")
  file.create(f_mine)
  db_register_file(f_mine, kind = "temp")

  f_theirs <- file.path(tempdir(), "someone_elses_file.tif")
  file.create(f_theirs)

  db_clean_temp(quiet = TRUE)
  expect_false(file.exists(f_mine))
  expect_true(file.exists(f_theirs))   # untouched
  unlink(f_theirs)
})

test_that("cleanup tracks provenance rather than inferring it from the path", {
  # An output written under tempdir() is still the user's file. Inferring
  # ownership from location would delete it.
  f <- file.path(tempdir(), "diamondback_output_in_tempdir.tif")
  file.create(f)
  db_register_file(f, kind = "output")
  db_clean_temp(quiet = TRUE)
  expect_true(file.exists(f))
  db_clean_temp(include_outputs = TRUE, quiet = TRUE)
  expect_false(file.exists(f))
})

test_that("db_clean_temp leaves user output files alone by default", {
  skip_if_no_python()
  f <- withr::local_tempfile(fileext = ".tif")
  label_patches(matrix(c(1, 0, 0, 1), 2, 2), output = f, quiet = TRUE)
  db_clean_temp(quiet = TRUE)
  expect_true(file.exists(f))   # an explicit output is the user's, not ours
})
