test_that("a three-step series produces one row per transition", {
  skip_if_no_python()
  m1 <- matrix(0, 6, 6); m1[2:3, 2:3] <- 1
  m2 <- matrix(0, 6, 6); m2[2:5, 2:5] <- 1
  m3 <- matrix(0, 6, 6); m3[2:3, 2:5] <- 1; m3[5, 2:5] <- 1

  s <- track_patch_series(list(m1, m2, m3), time = c(2000, 2010, 2020), quiet = TRUE)
  expect_s3_class(s, "patch_series")
  expect_equal(nrow(s$summary), 2L)
  expect_equal(s$summary$time_from, c(2000, 2010))
  expect_equal(s$summary$time_to, c(2010, 2020))
  expect_length(s$results, 3L)
  expect_named(s$results, c("2000", "2010", "2020"))
})

test_that("transitions carry time_from and time_to", {
  skip_if_no_python()
  m1 <- matrix(0, 5, 5); m1[2:4, 2:4] <- 1
  m2 <- matrix(0, 5, 5); m2[2, 2:4] <- 1; m2[4, 2:4] <- 1
  s <- track_patch_series(list(m1, m2), time = c(1985, 2024), quiet = TRUE)

  expect_true(all(c("time_from", "time_to") %in% names(s$transitions)))
  expect_true(all(s$transitions$time_from == 1985))
  expect_true(all(s$transitions$time_to == 2024))
  expect_true("split" %in% s$transitions$event)
  expect_equal(s$summary$split, 1L)
})

test_that("the series summary tallies events per transition", {
  skip_if_no_python()
  m1 <- matrix(0, 5, 5); m1[2:4, 2:4] <- 1
  m2 <- matrix(0, 5, 5); m2[2, 2:4] <- 1; m2[4, 2:4] <- 1   # split
  m3 <- matrix(0, 5, 5); m3[2:4, 2:4] <- 1                  # merge back
  s <- track_patch_series(list(m1, m2, m3), quiet = TRUE)

  expect_equal(s$summary$split, c(1L, 0L))
  expect_equal(s$summary$merger, c(0L, 1L))
  expect_equal(s$summary$n_patches_from, c(1L, 2L))
  expect_equal(s$summary$n_patches_to, c(2L, 1L))
})

test_that("a multi-layer SpatRaster works as a series", {
  skip_if_no_python()
  r1 <- rast_from(matrix(c(1, 1, 0, 0), 2, 2))
  r2 <- rast_from(matrix(c(1, 0, 0, 0), 2, 2))
  st <- c(r1, r2)
  names(st) <- c("y2000", "y2010")
  s <- track_patch_series(st, quiet = TRUE)
  expect_equal(nrow(s$summary), 1L)
  expect_equal(s$metadata$time, c("y2000", "y2010"))
})

test_that("filenames become default time labels", {
  skip_if_no_python()
  f1 <- tempfile("forest_1985_", fileext = ".tif")
  f2 <- tempfile("forest_2024_", fileext = ".tif")
  terra::writeRaster(rast_from(matrix(c(1, 1, 0, 0), 2, 2)), f1)
  terra::writeRaster(rast_from(matrix(c(1, 0, 0, 0), 2, 2)), f2)

  s <- track_patch_series(c(f1, f2), quiet = TRUE)
  expect_equal(s$metadata$time, tools::file_path_sans_ext(basename(c(f1, f2))))
})

test_that("a series with output_dir caches each step", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  f1 <- tempfile(fileext = ".tif"); f2 <- tempfile(fileext = ".tif")
  terra::writeRaster(rast_from(matrix(c(1, 1, 0, 0), 2, 2)), f1)
  terra::writeRaster(rast_from(matrix(c(1, 0, 0, 0), 2, 2)), f2)

  track_patch_series(c(f1, f2), time = c(2000, 2010), output_dir = d, quiet = TRUE)
  expect_true(file.exists(file.path(d, "2000", "patches.tif")))
  expect_true(file.exists(file.path(d, "2010", "metadata.json")))

  # Second run reuses both, which is what makes a long series resumable.
  expect_message(
    track_patch_series(c(f1, f2), time = c(2000, 2010), output_dir = d, quiet = FALSE),
    "Reusing cached result"
  )
})

test_that("a series of patch_result objects is accepted", {
  skip_if_no_python()
  p1 <- analyze_patches(matrix(c(1, 1, 0, 0), 2, 2), quiet = TRUE)
  p2 <- analyze_patches(matrix(c(1, 0, 0, 0), 2, 2), quiet = TRUE)
  s <- track_patch_series(list(p1, p2), time = c(1, 2), quiet = TRUE)
  expect_equal(nrow(s$summary), 1L)
})

test_that("a one-step series is rejected", {
  skip_if_no_python()
  expect_error(track_patch_series(list(matrix(1, 2, 2)), quiet = TRUE),
               "at least two")
  expect_error(track_patch_series(rast_from(matrix(1, 2, 2)), quiet = TRUE),
               "at least two")
})

test_that("mismatched time labels are rejected", {
  skip_if_no_python()
  expect_error(
    track_patch_series(list(matrix(1, 2, 2), matrix(1, 2, 2)), time = 1, quiet = TRUE),
    "one entry per time step"
  )
  expect_error(
    track_patch_series(list(matrix(1, 2, 2), matrix(1, 2, 2)), time = c(1, 1), quiet = TRUE),
    "duplicate"
  )
})

test_that("missing files in a series are reported before any work is done", {
  skip_if_no_python()
  expect_error(track_patch_series(c("no_such_a.tif", "no_such_b.tif"), quiet = TRUE),
               "do(es)? not exist")
})
