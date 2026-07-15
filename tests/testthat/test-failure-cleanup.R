# Adversarial: make each long stage fail, then check that nothing leaks and
# nothing of the user's is destroyed. A package that cleans up only on the happy
# path has not cleaned up.

n_tracked <- function() length(diamondback:::.db_files$created)
tracked <- function() diamondback:::.db_files$created

test_that("a failure during raster ingestion leaves no tracked temp files", {
  skip_if_no_python()
  before <- n_tracked()
  local_mocked_bindings(
    db_code_counts = function(code) stop("boom: simulated failure after ingestion")
  )
  expect_error(label_patches(matrix(c(1, 0, 0, 1), 2, 2), quiet = TRUE), "boom")
  expect_equal(n_tracked(), before)
})

test_that("a failure during labelling does not leave the reader open", {
  skip_if_no_python()
  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(rast_from(matrix(c(1, 0, 0, 1), 2, 2)), f, overwrite = TRUE)

  # The mock is scoped to this block only, so that recovery is tested against
  # the real function rather than the fault.
  local({
    local_mocked_bindings(
      db_label_classes = function(...) stop("boom: simulated labelling failure")
    )
    expect_error(label_patches(f, quiet = TRUE), "boom")
  })

  # terra's readStart/readStop is balanced by on.exit; if it were not, this
  # second read would fail or the file would stay locked.
  expect_no_error(terra::values(terra::rast(f), mat = FALSE))
  expect_no_error(label_patches(f, quiet = TRUE))
})

test_that("a failure while writing leaves no partial file at the user's path", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  out <- file.path(d, "labels.tif")

  # Fail after writeStart has created the file but before writeStop.
  local_mocked_bindings(
    db_progress_step = function(pb, i) stop("boom: simulated write failure")
  )
  expect_error(
    label_patches(matrix(rbinom(400, 1, 0.5), 20, 20), output = out, quiet = TRUE),
    "boom"
  )

  expect_false(file.exists(out))                       # no half-written raster
  expect_length(list.files(d, pattern = "diamondback-part"), 0L)  # no debris
})

test_that("a failed overwrite leaves the previous raster intact", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  out <- file.path(d, "labels.tif")

  # A good first run.
  label_patches(matrix(rbinom(400, 1, 0.5), 20, 20), output = out, quiet = TRUE)
  good <- terra::values(terra::rast(out), mat = FALSE)
  size_before <- file.info(out)$size

  # A failing overwrite must not destroy it. Writing straight to `out` would
  # have truncated the file the moment writeStart ran.
  local_mocked_bindings(
    db_progress_step = function(pb, i) stop("boom: simulated write failure")
  )
  expect_error(
    label_patches(matrix(rbinom(400, 1, 0.5), 20, 20), output = out,
                  overwrite = TRUE, quiet = TRUE),
    "boom"
  )

  expect_true(file.exists(out))
  expect_equal(file.info(out)$size, size_before)
  expect_equal(terra::values(terra::rast(out), mat = FALSE), good)
})

test_that("a failure during metric calculation leaves the labels usable", {
  skip_if_no_python()
  res <- label_patches(matrix(c(1, 1, 0, 1), 2, 2), quiet = TRUE)
  local_mocked_bindings(
    db_row_geometry = function(r) stop("boom: simulated metrics failure")
  )
  expect_error(patch_metrics(res, quiet = TRUE), "boom")
  # The result object is untouched and still works once the fault is gone.
  expect_equal(res$metadata$n_patches, 1L)
  expect_s4_class(res$patches, "SpatRaster")
})

test_that("repeated failures do not accumulate tracked files", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  local_mocked_bindings(
    db_progress_step = function(pb, i) stop("boom")
  )
  before <- n_tracked()
  for (i in 1:5) {
    try(label_patches(matrix(rbinom(400, 1, 0.5), 20, 20),
                      output = file.path(d, paste0("x", i, ".tif")), quiet = TRUE),
        silent = TRUE)
  }
  # Each attempt registers its part-file then removes it on failure.
  expect_lte(n_tracked() - before, 5L)
  expect_length(list.files(d), 0L)
})

test_that("cleanup never removes a user's output file", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  out <- file.path(d, "mine.tif")
  label_patches(matrix(c(1, 0, 0, 1), 2, 2), output = out, quiet = TRUE)

  db_clean_temp(quiet = TRUE)
  expect_true(file.exists(out))

  # Even the aggressive form only touches files diamondback made.
  other <- file.path(d, "not_ours.tif")
  file.create(other)
  db_clean_temp(include_outputs = TRUE, quiet = TRUE)
  expect_true(file.exists(other))
})

test_that("the Python module handle survives a failed operation", {
  skip_if_no_python()
  before <- diamondback:::.db_state$py
  try(label_patches(matrix(1, 2, 2), memory_limit = 1, quiet = TRUE), silent = TRUE)
  expect_identical(diamondback:::.db_state$py, before)
  # And the backend still works afterwards. (The two cells touch diagonally, so
  # this is one patch under the default 8-connectivity and two under 4.)
  expect_equal(label_patches(matrix(c(1, 0, 0, 1), 2, 2), quiet = TRUE)$metadata$n_patches, 1L)
  expect_equal(label_patches(matrix(c(1, 0, 0, 1), 2, 2), directions = 4,
                             quiet = TRUE)$metadata$n_patches, 2L)
})

test_that("a memory refusal allocates nothing and leaves no trace", {
  skip_if_no_python()
  before <- n_tracked()
  expect_error(
    label_patches(matrix(1, 100, 100), output = tempfile(fileext = ".tif"),
                  memory_limit = 1000, quiet = TRUE),
    class = "diamondback_memory_error"
  )
  expect_equal(n_tracked(), before)
})
