# Geometry the package does not support must be rejected, not approximated.

test_that("a lon/lat raster reaching beyond the pole is rejected", {
  g <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 10,
                   ymin = 80, ymax = 100, crs = "EPSG:4326")
  terra::values(g) <- 1
  # terra builds this happily; the geodesic row geometry would be nonsense.
  expect_error(label_patches(g, quiet = TRUE), "beyond the poles")
})

test_that("a lon/lat raster at exactly the poles is allowed", {
  g <- terra::rast(nrows = 4, ncols = 4, xmin = -180, xmax = 180,
                   ymin = -90, ymax = 90, crs = "EPSG:4326")
  terra::values(g) <- 1
  expect_no_error(label_patches(g, quiet = TRUE))
})

test_that("a projected raster with a large extent is not mistaken for lat/lon", {
  r <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4e6,
                   ymin = 0, ymax = 4e6, crs = "EPSG:5070")
  terra::values(r) <- 1
  expect_no_error(label_patches(r, quiet = TRUE))
})

test_that("rotated rasters are rejected with a route forward", {
  # terra's data model has no rotation term, so this cannot be built directly;
  # mock is.rotated to pin the guard itself.
  r <- rast_from(matrix(1, 3, 3))
  local_mocked_bindings(db_is_lonlat = function(r) FALSE)
  testthat::local_mocked_bindings(is.rotated = function(x) TRUE, .package = "terra")
  expect_error(db_check_grid(r), "rotated raster")
  expect_error(db_check_grid(r), "rectify")
})

test_that("a normal raster passes the grid check", {
  expect_true(db_check_grid(rast_from(matrix(1, 3, 3))))
  g <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3,
                   crs = "EPSG:4326")
  terra::values(g) <- 1
  expect_true(db_check_grid(g))
})

# --- class encoding ---------------------------------------------------------

test_that("more than 253 classes works via a wider code array", {
  skip_if_no_python()
  # 300 distinct values: the old uint8 encoding (3 + k) could not hold these.
  n <- 300
  m <- matrix(seq_len(n), nrow = 20, ncol = 15)
  r <- rast_from(m, res = 1)

  # as.list(): 300 separate classes. seq_len(n) would be one class of 300 values.
  res <- label_patches(r, class = as.list(seq_len(n)), quiet = TRUE)
  expect_equal(res$metadata$n_patches, n)
  expect_equal(sort(as.numeric(unique(res$metrics$class))), seq_len(n))
  expect_equal(sum(res$metrics$cells), n)
  expect_true(validate_patch_result(res, quiet = TRUE)$ok)
})

test_that("the code dtype widens exactly at the uint8 boundary", {
  skip_if_no_python()
  py <- diamondback:::db_py()
  expect_equal(py_r(py$code_dtype(1L)$name), "uint8")
  expect_equal(py_r(py$code_dtype(253L)$name), "uint8")
  expect_equal(py_r(py$code_dtype(254L)$name), "uint16")
  expect_equal(py_r(py$code_dtype(65533L)$name), "uint16")
  expect_error(py_r(py$code_dtype(65534L)), "at most")
})

test_that("the memory estimate accounts for the wider code array", {
  n <- 1e6
  expect_equal(db_estimate_memory(n, "label", code_bytes = 1), n * 6)
  expect_equal(db_estimate_memory(n, "label", code_bytes = 2), n * 7)
})

test_that("too many classes is rejected with a clear limit", {
  skip_if_no_python()
  expect_error(
    label_patches(matrix(1, 2, 2), class = as.list(seq_len(70000)), quiet = TRUE),
    "at most"
  )
})

test_that("class = 'all' warns when the value count implies many passes", {
  skip_if_no_python()
  m <- matrix(seq_len(150), nrow = 10, ncol = 15)
  expect_warning(label_patches(rast_from(m, res = 1), class = "all", quiet = TRUE),
                 "separate pass")
})

test_that("class = 'all' rejects a continuous raster politely", {
  skip_if_no_python()
  set.seed(9)
  r <- rast_from(matrix(runif(70000), nrow = 250, ncol = 280), res = 1)
  expect_error(label_patches(r, class = "all", quiet = TRUE), "continuous rather than categorical")
})

# --- fingerprint modes ------------------------------------------------------

test_that("fingerprint = 'full' hashes a file that 'auto' would skip", {
  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(rast_from(matrix(1, 4, 4)), f, overwrite = TRUE)

  auto <- db_source_info(f, "auto")
  full <- db_source_info(f, "full")
  fast <- db_source_info(f, "fast")

  expect_false(is.na(auto$hash))          # small file: auto hashes it
  expect_equal(auto$hash, full$hash)
  expect_true(is.na(fast$hash))           # fast never hashes
  expect_equal(fast$fingerprint, "fast")
})

test_that("'auto' falls back to size and mtime for a large file", {
  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(rast_from(matrix(1, 4, 4)), f, overwrite = TRUE)
  # Pretend the file is huge without writing 300 MB.
  info <- db_source_info(f, "auto", hash_max_mb = 0)
  expect_true(is.na(info$hash))
  expect_false(is.na(info$size))
  expect_false(is.na(info$mtime))
})

test_that("the fingerprint mode is part of the cache key", {
  skip_if_no_python()
  f <- withr::local_tempfile(fileext = ".tif")
  d <- withr::local_tempdir()
  terra::writeRaster(rast_from(matrix(c(1, 0, 0, 1), 2, 2)), f, overwrite = TRUE)

  analyze_patches(f, output_dir = d, fingerprint = "fast", quiet = TRUE)

  # A result cached on size+mtime must not satisfy a request asking for a
  # content hash, and the reported reason must be the mode, not the hash it
  # caused to differ.
  msgs <- capture.output(
    analyze_patches(f, output_dir = d, fingerprint = "full", quiet = FALSE),
    type = "message")
  expect_match(paste(msgs, collapse = " "), "source\\$fingerprint differs")

  # The same mode then hits.
  msgs2 <- capture.output(
    analyze_patches(f, output_dir = d, fingerprint = "full", quiet = FALSE),
    type = "message")
  expect_match(paste(msgs2, collapse = " "), "Reusing cached result")
})

test_that("a content change with identical size and mtime is caught by the hash", {
  # The claim under test is about the cache key, not about the filesystem, so it
  # is tested on the key directly. Doing it by touching timestamps depends on
  # Sys.setFileTime() actually restoring an mtime exactly, which is not true
  # everywhere -- on Windows CI the mtime survived only approximately, the
  # comparison reported "mtime differs", and the test failed for a reason that
  # had nothing to do with the behaviour it was checking.
  base <- list(type = "file", path = "/x/y.tif", size = 1234, mtime = "2026-01-01 00:00:00",
               hash = "aaaa", fingerprint = "full")
  meta <- function(src) {
    list(algorithm_version = "1", source = src, mask_source = NULL,
         geometry = list(nrow = 2, ncol = 2, xmin = 0, xmax = 20, ymin = 0,
                         ymax = 20, xres = 10, yres = 10, crs = "x"),
         class = NA_real_, binary = TRUE, directions = 8, na = "outside",
         crop = FALSE)
  }

  same <- meta(base)
  expect_true(isTRUE(db_cache_match(same, same)))

  # Identical size and mtime, different contents: only the hash can tell.
  changed <- base; changed$hash <- "bbbb"
  expect_equal(db_cache_match(meta(base), meta(changed)), "source$hash")

  # And under "fast" there is no hash to compare, so such a change is invisible.
  # That is the documented trade-off, pinned here so it cannot drift silently.
  fast_a <- base; fast_a$fingerprint <- "fast"; fast_a$hash <- NA_character_
  fast_b <- fast_a
  expect_true(isTRUE(db_cache_match(meta(fast_a), meta(fast_b))))
})

test_that("fingerprint = 'full' rereads contents when a file is rewritten", {
  skip_if_no_python()
  f <- withr::local_tempfile(fileext = ".tif")
  d <- withr::local_tempdir()

  r <- rast_from(matrix(c(1, 0, 0, 1), 2, 2))
  terra::writeRaster(r, f, overwrite = TRUE)
  analyze_patches(f, output_dir = d, fingerprint = "full", quiet = TRUE)
  mt <- file.info(f)$mtime

  # Same dimensions and dtype -> same file size; try to restore the timestamp.
  r2 <- rast_from(matrix(c(1, 1, 1, 1), 2, 2))
  terra::writeRaster(r2, f, overwrite = TRUE)
  Sys.setFileTime(f, mt)

  # Not every filesystem restores an mtime exactly. When it does, the hash is
  # the only thing that can notice the change, and that is what we check. When
  # it does not, the premise is gone and the test would be checking nothing.
  restored <- isTRUE(all.equal(as.character(file.info(f)$mtime), as.character(mt)))
  skip_if_not(restored, "filesystem did not restore the mtime exactly")

  msgs <- capture.output(
    analyze_patches(f, output_dir = d, fingerprint = "full", quiet = FALSE),
    type = "message")
  expect_match(paste(msgs, collapse = " "), "source\\$hash differs")
})
