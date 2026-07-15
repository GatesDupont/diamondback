# Patch IDs have to survive SciPy -> NumPy -> reticulate -> R -> GDAL -> back.
# The narrowest link in that chain, not SciPy, sets the ceiling.

test_that("terra preserves extreme patch IDs exactly through a write and reopen", {
  # This is the round trip every result takes. INT4S plus INT_MIN as the NA flag
  # must reproduce 1, a large ID, the maximum, and NA without drift.
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 20, ymin = 0, ymax = 20,
                   crs = "EPSG:5070")
  vals <- c(1L, 123456789L, .Machine$integer.max, NA_integer_)
  terra::values(r) <- vals

  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE, datatype = "INT4S",
                     NAflag = -2147483648)
  back <- terra::values(terra::rast(f), mat = FALSE)

  expect_identical(as.integer(back), vals)
  expect_true(is.na(back[4]))
  expect_equal(back[3], .Machine$integer.max)
})

test_that("the NA sentinel is not confusable with a real patch ID", {
  # INT_MIN is R's NA_integer_ and the raster's no-data value, so the usable ID
  # range must stop one short of it in magnitude.
  expect_equal(MAX_PATCH_ID, .Machine$integer.max)
  expect_true(MAX_PATCH_ID < abs(-2147483648))
})

test_that("a patch count beyond the representable range is refused, not overflowed", {
  expect_no_error(db_check_label_range(MAX_PATCH_ID))
  expect_error(db_check_label_range(MAX_PATCH_ID + 1),
               class = "diamondback_overflow_error")
  err <- tryCatch(db_check_label_range(3e9), error = function(e) conditionMessage(e))
  expect_match(err, "more than can be represented")
  expect_match(err, "2,147,483,647")
})

test_that("labels round-trip exactly through the real write path", {
  skip_if_no_python()
  # Not extreme values, but the actual pipeline: every label written must read
  # back identical, and NA must stay NA rather than becoming 0 or a huge number.
  set.seed(3)
  m <- matrix(rbinom(400, 1, 0.4), 20, 20)
  m[sample(400, 40)] <- NA

  f <- withr::local_tempfile(fileext = ".tif")
  disk <- label_patches(rast_from(m), output = f, quiet = TRUE)
  mem <- label_patches(rast_from(m), quiet = TRUE)

  dv <- terra::values(terra::rast(f), mat = FALSE)
  mv <- terra::values(mem$patches, mat = FALSE)
  expect_identical(dv, mv)
  expect_equal(sum(!is.na(dv) & dv > 0), disk$metadata$cells$foreground)
  expect_equal(sum(is.na(dv)), disk$metadata$cells$missing + disk$metadata$cells$outside)
  expect_equal(max(dv, na.rm = TRUE), disk$metadata$n_patches)
})

test_that("the label array dtype widens with cell count, not patch count", {
  skip_if_no_python()
  # int32 is chosen while the cell count fits it: a patch needs at least one
  # cell, so patch IDs cannot exceed the cell count.
  py <- diamondback:::db_py()
  code <- py$code_block(as.numeric(c(1, 0, 0, 1)), 2L, 2L, mask = NULL,
                        class_values = NULL, n_classes = 1L)

  small <- py$label_array(code, 0L, 8L, FALSE)
  expect_equal(py_r(py_get(small, "labels")$dtype$name), "int32")

  big <- py$label_array(code, 0L, 8L, TRUE)
  expect_equal(py_r(py_get(big, "labels")$dtype$name), "int64")
})

test_that("an int64 label array still returns int32 cells to R", {
  skip_if_no_python()
  # A raster over 2^31 cells labels into int64; the values handed back to R must
  # still be int32, because that is what an R integer and a GeoTIFF band hold.
  skip_if_not(diamondback_ready())
  py <- diamondback:::db_py()
  code <- py$code_block(as.numeric(c(1, 0, 0, 1)), 2L, 2L, mask = NULL,
                        class_values = NULL, n_classes = 1L)
  lab <- py_get(py$label_array(code, 0L, 8L, TRUE), "labels")
  rows <- py$output_rows(lab, code, 0L, 2L, FALSE)
  expect_equal(py_r(rows$dtype$name), "int32")
  expect_type(as.integer(py_r(rows)), "integer")
})

test_that("validation reports the ID range as a check", {
  skip_if_no_python()
  v <- validate_patch_result(analyze_patches(matrix(c(1, 0, 0, 1), 2, 2), quiet = TRUE),
                             quiet = TRUE)
  expect_true(any(grepl("32-bit", v$checks$check)))
  expect_true(all(v$checks$passed))
})

test_that("overlap cross-tabulation refuses counts that would overflow its key", {
  skip_if_no_python()
  # The crosstab packs (id1, id2) into one int64 key. Rather than wrap around
  # silently it raises, and that surfaces as a diamondback error.
  py <- diamondback:::db_py()
  code <- py$code_block(as.numeric(c(1, 0, 0, 1)), 2L, 2L, mask = NULL,
                        class_values = NULL, n_classes = 1L)
  lab <- py_get(py$label_array(code, 0L, 8L, FALSE), "labels")
  expect_error(
    py_try(py$overlap_counts(lab, lab, 4e9, 4e9), "cross-tabulating"),
    class = "diamondback_python_error"
  )
})
