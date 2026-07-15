# One element of `class` is one class. A numeric vector is a single class made
# of several raster values; a list is several classes. This is what removes the
# need to hand-build a binary raster first -- the step where NA gets destroyed.

test_that("a numeric vector is ONE class, connected across its values", {
  skip_if_no_python()
  # 41 and 42 are adjacent. As one class they are one patch; the boundary
  # between two values inside a group is not an edge.
  m <- matrix(c(41, 42, 21, 43), nrow = 2)
  res <- label_patches(rast_from(m), class = c(41, 42, 43), quiet = TRUE)
  expect_equal(res$metadata$n_patches, 1L)
  expect_equal(res$metrics$cells, 3)
  expect_false("class" %in% names(res$metrics))   # one class needs no column
})

test_that("a list is several classes, never joined to each other", {
  skip_if_no_python()
  m <- matrix(c(41, 42, 21, 43), nrow = 2)
  res <- label_patches(rast_from(m), class = list(41, 42, 43), quiet = TRUE)
  expect_equal(res$metadata$n_patches, 3L)
  expect_equal(sort(res$metrics$class), c("41", "42", "43"))
})

test_that("adjacency within a group is habitat; adjacency between groups is not", {
  skip_if_no_python()
  # A solid 2x2 of forest values touching a solid 2x2 of wetland values.
  m <- matrix(c(41, 42, 90, 95,
                42, 41, 95, 90), nrow = 2, byrow = TRUE)
  res <- label_patches(rast_from(m),
                       class = list(forest = c(41, 42), wetland = c(90, 95)),
                       directions = 8, quiet = TRUE)
  # Each group is internally continuous -> one patch each, not four.
  expect_equal(res$metadata$n_patches, 2L)
  expect_equal(sort(res$metrics$cells), c(4, 4))
  expect_setequal(res$metrics$class, c("forest", "wetland"))

  # And the two classes are not merged despite touching.
  expect_equal(nrow(res$metrics), 2L)
})

test_that("named groups reach metrics and metadata", {
  skip_if_no_python()
  m <- matrix(c(41, 42, 90, 95), nrow = 2)
  res <- label_patches(rast_from(m),
                       class = list(forest = c(41, 42), wetland = c(90, 95)),
                       quiet = TRUE)
  expect_setequal(res$metrics$class, c("forest", "wetland"))
  expect_setequal(res$metadata$class$labels, c("forest", "wetland"))
  expect_equal(sort(unname(unlist(res$metadata$class$groups))), c(41, 42, 90, 95))
  expect_false(isTRUE(res$metadata$binary))
})

test_that("unnamed groups get stable labels from their values", {
  skip_if_no_python()
  m <- matrix(c(41, 42, 90, 95), nrow = 2)
  res <- label_patches(rast_from(m), class = list(c(41, 42), c(90, 95)), quiet = TRUE)
  # "41+42" says what the class is; "class_1" would not.
  expect_setequal(res$metrics$class, c("41+42", "90+95"))
})

test_that("partially named groups fill the gaps from values", {
  skip_if_no_python()
  m <- matrix(c(41, 42, 90, 95), nrow = 2)
  res <- label_patches(rast_from(m), class = list(forest = c(41, 42), c(90, 95)),
                       quiet = TRUE)
  expect_setequal(res$metrics$class, c("forest", "90+95"))
})

test_that("duplicate values inside a group are harmless", {
  skip_if_no_python()
  m <- matrix(c(41, 42, 21, 43), nrow = 2)
  a <- label_patches(rast_from(m), class = c(41, 41, 42, 43, 43), quiet = TRUE)
  b <- label_patches(rast_from(m), class = c(41, 42, 43), quiet = TRUE)
  expect_equal(a$metadata$n_patches, b$metadata$n_patches)
  expect_equal(a$metrics$cells, b$metrics$cells)
  expect_equal(unname(unlist(a$metadata$class$groups)), c(41, 42, 43))  # normalised
})

test_that("a value in two groups is rejected", {
  skip_if_no_python()
  # A cell carries one code, so this has no consistent answer.
  expect_error(
    label_patches(matrix(41, 2, 2), class = list(a = c(41, 42), b = c(42, 43)), quiet = TRUE),
    "more than one class"
  )
  err <- tryCatch(
    label_patches(matrix(41, 2, 2), class = list(c(41, 42), c(42, 43)), quiet = TRUE),
    error = conditionMessage)
  expect_match(err, "42")     # names the offending value
})

test_that("empty and NA groups are rejected", {
  skip_if_no_python()
  expect_error(label_patches(matrix(1, 2, 2), class = list(), quiet = TRUE), "empty")
  expect_error(label_patches(matrix(1, 2, 2), class = list(a = numeric(0)), quiet = TRUE),
               "is empty")
  expect_error(label_patches(matrix(1, 2, 2), class = list(a = c(1, NA)), quiet = TRUE), "NA")
  expect_error(label_patches(matrix(1, 2, 2), class = list(a = NA), quiet = TRUE), "NA")
  expect_error(label_patches(matrix(1, 2, 2), class = c(1, NA), quiet = TRUE), "NA")
  expect_error(label_patches(matrix(1, 2, 2), class = NA, quiet = TRUE), "NA")
})

test_that("duplicate class names are rejected", {
  skip_if_no_python()
  expect_error(
    label_patches(matrix(1, 2, 2), class = list(a = 1, a = 2), quiet = TRUE),
    "duplicate class name"
  )
})

test_that("non-numeric groups are rejected clearly", {
  skip_if_no_python()
  expect_error(label_patches(matrix(1, 2, 2), class = list(a = "x"), quiet = TRUE),
               "must be numeric")
})

test_that("class = 'all' still gives one class per value", {
  skip_if_no_python()
  m <- matrix(c(1, 1, 2, 3), nrow = 2)
  res <- label_patches(rast_from(m), class = "all", directions = 4, quiet = TRUE)
  expect_setequal(res$metrics$class, c("1", "2", "3"))
  expect_equal(sum(res$metrics$cells), 4)
})

test_that("grouping removes the need for a lossy pre-reclass", {
  skip_if_no_python()
  # The bug this feature exists to prevent. The hand-rolled route --
  # ifel(x %in% classes, 1, 0) -- sends NA to 0 and silently converts unknown
  # data into valid background. Grouping never sees a reclassified raster.
  m <- matrix(c(41, 42, NA, 21, 43, NA, 21, 21, 41), nrow = 3, byrow = TRUE)
  r <- rast_from(m)

  grouped <- label_patches(r, class = c(41, 42, 43), quiet = TRUE)
  expect_equal(grouped$metadata$cells$missing, 2)      # NA kept as its own state
  expect_equal(grouped$metadata$cells$foreground, 4)

  # The hand-rolled reclass, done on the values so the point is explicit:
  # `%in%` sends NA to FALSE, which ifel() then writes as 0.
  vals <- terra::values(r, mat = FALSE)
  flat <- r
  terra::values(flat) <- as.integer(vals %in% c(41, 42, 43))   # NA -> FALSE -> 0
  lossy <- label_patches(flat, class = 1, quiet = TRUE)
  expect_equal(lossy$metadata$cells$missing, 0)        # the NA is gone
  expect_gt(lossy$metadata$cells$background, grouped$metadata$cells$background)
})

test_that("the cache tells one grouped class apart from several separate ones", {
  skip_if_no_python()
  f <- withr::local_tempfile(fileext = ".tif")
  d <- withr::local_tempdir()
  terra::writeRaster(rast_from(matrix(c(41, 42, 21, 43), nrow = 2)), f, overwrite = TRUE)

  analyze_patches(f, class = c(41, 42, 43), output_dir = d, quiet = TRUE)
  # Same values, different meaning: this must not be served the grouped result.
  expect_message(analyze_patches(f, class = list(41, 42, 43), output_dir = d, quiet = FALSE),
                 "class differs")
  expect_message(analyze_patches(f, class = list(41, 42, 43), output_dir = d, quiet = FALSE),
                 "Reusing cached result")
})

test_that("group order does not change the cache key", {
  # The key is sorted by label, so writing the same classes in another order is
  # the same computation and should hit.
  mk <- function(cl) list(algorithm_version = "2", binary = FALSE, class = cl,
                          source = list(type = "file", path = "p", size = 1,
                                        mtime = "t", hash = "h", fingerprint = "auto"),
                          mask_source = NULL,
                          geometry = list(nrow = 1, ncol = 1, xmin = 0, xmax = 1,
                                          ymin = 0, ymax = 1, xres = 1, yres = 1, crs = ""),
                          directions = 8, na = "outside", crop = FALSE)
  a <- mk(list(labels = c("forest", "wet"), groups = list(c(41, 42), c(90, 95))))
  b <- mk(list(labels = c("wet", "forest"), groups = list(c(90, 95), c(41, 42))))
  expect_true(isTRUE(db_cache_match(a, b)))

  # But a genuinely different grouping must not collide.
  cc <- mk(list(labels = c("41+42+90+95"), groups = list(c(41, 42, 90, 95))))
  expect_false(isTRUE(db_cache_match(a, cc)))
})
