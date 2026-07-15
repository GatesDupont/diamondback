test_that("a healthy result passes every check", {
  skip_if_no_python()
  set.seed(41)
  m <- matrix(rbinom(400, 1, 0.4), 20, 20)
  v <- validate_patch_result(analyze_patches(rast_from(m), quiet = TRUE), quiet = TRUE)
  expect_true(v$ok)
  expect_true(all(v$checks$passed))
  expect_gt(nrow(v$checks), 5)
})

test_that("validation passes with NA, masks and multiple classes", {
  skip_if_no_python()
  v <- matrix(c(1, 2, NA, 1, 1, 2, NA, 2, 1), nrow = 3)
  mk <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 0), nrow = 3)
  res <- analyze_patches(rast_from(v), class = c(1, 2), mask = rast_from(mk),
                         crop = FALSE, quiet = TRUE)
  expect_true(validate_patch_result(res, quiet = TRUE)$ok)
})

test_that("validation passes on an empty result", {
  skip_if_no_python()
  expect_warning(res <- analyze_patches(matrix(0, 4, 4), quiet = TRUE))
  expect_true(validate_patch_result(res, quiet = TRUE)$ok)
})

test_that("validation catches a metrics table that disagrees with the raster", {
  skip_if_no_python()
  m <- matrix(0, 6, 6); m[1, 1] <- 1; m[3, 3] <- 1; m[5, 5] <- 1
  res <- analyze_patches(rast_from(m), directions = 4, quiet = TRUE)

  # Corrupt the table the way a bad join would.
  bad <- res
  bad$metrics <- bad$metrics[1:2, ]
  v <- validate_patch_result(bad, quiet = TRUE)
  expect_false(v$ok)
  expect_true(any(grepl("metrics rows", v$checks$check[!v$checks$passed])))
})

test_that("validation catches impossible core area", {
  skip_if_no_python()
  m <- matrix(0, 6, 6); m[2:5, 2:5] <- 1
  res <- patch_core_area(analyze_patches(rast_from(m, res = 1), quiet = TRUE),
                         depth = 1, quiet = TRUE)
  bad <- res
  bad$metrics$core_cells <- bad$metrics$cells + 1
  v <- validate_patch_result(bad, quiet = TRUE)
  expect_false(v$ok)
  expect_true(any(grepl("core cells", v$checks$check[!v$checks$passed])))
})

test_that("validation catches an inconsistent edge decomposition", {
  skip_if_no_python()
  m <- matrix(0, 6, 6); m[2:4, 2:4] <- 1
  res <- analyze_patches(rast_from(m), quiet = TRUE)
  bad <- res
  bad$metrics$edge_valid <- bad$metrics$edge_valid + 1
  v <- validate_patch_result(bad, quiet = TRUE)
  expect_false(v$ok)
  expect_true(any(grepl("edge components", v$checks$check[!v$checks$passed])))
})

test_that("error = TRUE turns a failure into a condition", {
  skip_if_no_python()
  res <- analyze_patches(matrix(c(1, 0, 0, 1), 2, 2), quiet = TRUE)
  bad <- res
  bad$metrics <- bad$metrics[0, ]
  expect_error(validate_patch_result(bad, error = TRUE, quiet = TRUE),
               class = "diamondback_validation_error")
})

test_that("label_patches(validate = TRUE) records the outcome", {
  skip_if_no_python()
  res <- label_patches(matrix(c(1, 0, 0, 1), 2, 2), validate = TRUE, quiet = TRUE)
  expect_true(res$metadata$validated)
})

test_that("foreground never leaks outside the mask", {
  skip_if_no_python()
  set.seed(51)
  v <- matrix(rbinom(400, 1, 0.6), 20, 20)
  mk <- matrix(rbinom(400, 1, 0.7), 20, 20)
  res <- analyze_patches(rast_from(v), mask = rast_from(mk), crop = FALSE, quiet = TRUE)

  v_check <- validate_patch_result(res, quiet = TRUE)
  expect_true(v_check$ok)

  # Independently: every labelled cell must be inside the mask.
  lab <- terra::values(res$patches, mat = FALSE)
  mkv <- as.vector(t(mk))
  expect_true(all(is.na(lab[mkv == 0])))
  expect_equal(sum(!is.na(lab) & lab > 0), res$metadata$cells$foreground)
})
