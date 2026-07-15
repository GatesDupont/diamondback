test_that("$patches works whether the raster is in memory or on disk", {
  skip_if_no_python()
  m <- matrix(c(1, 1, 0, 0), 2, 2)

  mem <- label_patches(m, quiet = TRUE)
  expect_s4_class(mem$patches, "SpatRaster")
  expect_null(mem$path)

  f <- withr::local_tempfile(fileext = ".tif")
  disk <- label_patches(m, output = f, quiet = TRUE)
  expect_s4_class(disk$patches, "SpatRaster")
  expect_equal(disk$path, f)

  # Callers cannot tell the difference, which is the point.
  expect_equal(terra::values(mem$patches, mat = FALSE),
               terra::values(disk$patches, mat = FALSE))
})

test_that("[[ and $ agree on patches", {
  skip_if_no_python()
  res <- label_patches(matrix(c(1, 0, 0, 1), 2, 2), quiet = TRUE)
  expect_equal(terra::values(res$patches, mat = FALSE),
               terra::values(res[["patches"]], mat = FALSE))
})

test_that("other list elements are reachable normally", {
  skip_if_no_python()
  res <- analyze_patches(matrix(c(1, 0, 0, 1), 2, 2), quiet = TRUE)
  expect_s3_class(res$metrics, "data.frame")
  expect_type(res$metadata, "list")
  expect_true(is.null(res$nonexistent))
})

test_that("a missing label file gives a clear error rather than a terra crash", {
  skip_if_no_python()
  f <- withr::local_tempfile(fileext = ".tif")
  res <- label_patches(matrix(c(1, 0, 0, 1), 2, 2), output = f, quiet = TRUE)
  unlink(f)
  expect_error(res$patches, "no longer|missing")
})

# cli writes to stderr, so these are messages rather than stdout output.
shown <- function(expr) {
  paste(utils::capture.output(expr, type = "message"), collapse = "\n")
}

test_that("print and summary report the essentials", {
  skip_if_no_python()
  set.seed(61)
  m <- matrix(rbinom(400, 1, 0.4), 20, 20)
  res <- analyze_patches(rast_from(m), quiet = TRUE)

  out <- shown(print(res))
  expect_match(out, "patch_result")
  expect_match(out, "background")
  expect_match(out, "outside domain")   # the four states are always visible

  expect_match(shown(summary(res)), "Patch size")

  cmp <- compare_patches(res, res, quiet = TRUE)
  expect_match(shown(print(cmp)), "patch_comparison")
  expect_match(shown(summary(cmp)), "Area accounting")
})

test_that("print works on an empty result", {
  skip_if_no_python()
  expect_warning(res <- analyze_patches(matrix(0, 3, 3), quiet = TRUE))
  expect_match(shown(print(res)), "0 patches")
  expect_match(shown(summary(res)), "No patches")
})

test_that("patch_domain reconstructs all four cell states", {
  skip_if_no_python()
  v <- matrix(c(1, 0, NA,
                1, 0, 0,
                0, 0, 1), nrow = 3, byrow = TRUE)
  mk <- matrix(1, 3, 3); mk[3, 3] <- 0
  res <- label_patches(rast_from(v), mask = rast_from(mk), crop = FALSE, quiet = TRUE)

  dom <- patch_domain(res)
  expect_s4_class(dom, "SpatRaster")
  vals <- terra::values(dom, mat = FALSE)

  expect_equal(as.integer(vals[3]), 1L)   # NA in source -> missing
  expect_equal(as.integer(vals[9]), 0L)   # masked out -> outside
  expect_equal(as.integer(vals[1]), 3L)   # foreground -> patch
  expect_equal(as.integer(vals[2]), 2L)   # background

  # And the counts agree with the metadata, which is the cheap path to the
  # same facts.
  cl <- res$metadata$cells
  expect_equal(sum(vals == 0), cl$outside)
  expect_equal(sum(vals == 1), cl$missing)
  expect_equal(sum(vals == 2), cl$background)
  expect_equal(sum(vals == 3), cl$foreground)
})

test_that("analyze_patches equals label_patches plus patch_metrics", {
  skip_if_no_python()
  set.seed(71)
  m <- rast_from(matrix(rbinom(400, 1, 0.4), 20, 20))
  a <- analyze_patches(m, quiet = TRUE)
  b <- patch_metrics(label_patches(m, quiet = TRUE), quiet = TRUE)
  expect_equal(a$metrics, b$metrics)
})

test_that("analyze_patches(metrics = FALSE) returns counts only", {
  skip_if_no_python()
  res <- analyze_patches(matrix(c(1, 0, 0, 1), 2, 2), metrics = FALSE, quiet = TRUE)
  expect_named(res$metrics, c("patch_id", "cells"))
})

test_that("edge_depth in patch_metrics adds core columns", {
  skip_if_no_python()
  m <- matrix(0, 9, 9); m[2:8, 2:8] <- 1
  res <- patch_metrics(label_patches(rast_from(m, res = 1), quiet = TRUE),
                       edge_depth = 1, quiet = TRUE)
  expect_true(all(c("core_cells", "core_area_m2", "core_fraction") %in% names(res$metrics)))
  expect_equal(res$metrics$core_cells, 25)
})
