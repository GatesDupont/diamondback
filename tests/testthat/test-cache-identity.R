# What does a cache key mean for a SpatRaster that only exists in memory?
# There is no path, size or mtime to lean on, so the answer has to be either
# "hash the cells" or "decline to cache". It must never be "assume unchanged".

test_that("terra drops the file source once a raster is edited in memory", {
  # The whole in-memory story rests on this terra behaviour, so it is pinned
  # here: if terra ever kept the source after an edit, db_source_info() would
  # fingerprint the file while the data in hand differed from it, which is a
  # false cache hit waiting to happen.
  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(rast_from(matrix(c(1, 0, 0, 1), 2, 2)), f, overwrite = TRUE)

  r <- terra::rast(f)
  expect_equal(db_source_info(r)$type, "file")

  r[1, 1] <- 99
  expect_equal(db_source_info(r)$type, "memory")
  expect_false(is.na(db_source_info(r)$hash))
})

test_that("an edited in-memory raster gets a different fingerprint", {
  f <- withr::local_tempfile(fileext = ".tif")
  terra::writeRaster(rast_from(matrix(c(1, 0, 0, 1), 2, 2)), f, overwrite = TRUE)

  r1 <- terra::rast(f); terra::values(r1) <- c(1, 0, 0, 1)
  r2 <- terra::rast(f); terra::values(r2) <- c(1, 1, 0, 1)
  expect_false(identical(db_source_info(r1)$hash, db_source_info(r2)$hash))

  # Identical values give identical fingerprints, so caching still works.
  r3 <- terra::rast(f); terra::values(r3) <- c(1, 0, 0, 1)
  expect_identical(db_source_info(r1)$hash, db_source_info(r3)$hash)
})

test_that("derived rasters are fingerprinted by their own values, not the parent's", {
  f <- withr::local_tempfile(fileext = ".tif")
  r0 <- rast_from(matrix(1:16, 4, 4))
  terra::writeRaster(r0, f, overwrite = TRUE)
  r <- terra::rast(f)

  parent <- db_source_info(r)$hash
  for (derived in list(r * 2, terra::ifel(r > 8, 1, 0), terra::subst(r, 1, 99))) {
    info <- db_source_info(derived)
    expect_equal(info$type, "memory")
    expect_false(identical(info$hash, parent))
  }
})

test_that("an in-memory raster edited between runs misses the cache", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  r <- rast_from(matrix(c(1, 0, 0, 1), 2, 2))

  analyze_patches(r, output_dir = d, quiet = TRUE)
  expect_message(analyze_patches(r, output_dir = d, quiet = FALSE),
                 "Reusing cached result")

  # Change one cell in memory; the cache must not be reused.
  r[1, 1] <- 0
  expect_message(analyze_patches(r, output_dir = d, quiet = FALSE),
                 "source\\$hash differs")
})

test_that("fingerprint modes mean the right thing for in-memory rasters", {
  r <- rast_from(matrix(c(1, 0, 0, 1), 2, 2))
  # "full" always hashes; "auto" hashes while cheap; "fast" never does, so it
  # leaves an in-memory raster with no identity at all.
  expect_false(is.na(db_source_info(r, "full")$hash))
  expect_false(is.na(db_source_info(r, "auto")$hash))
  expect_true(is.na(db_source_info(r, "fast")$hash))
})

test_that("a raster too large to hash under 'auto' is never served from cache", {
  # No hash means no identity: matching on nothing would be a guess.
  fake <- list(type = "memory", path = NA_character_, size = 5e6,
               mtime = NA_character_, hash = NA_character_, fingerprint = "auto")
  meta <- list(algorithm_version = "1", source = fake, mask_source = NULL,
               geometry = list(nrow = 1, ncol = 1, xmin = 0, xmax = 1, ymin = 0,
                               ymax = 1, xres = 1, yres = 1, crs = ""),
               class = NA_real_, binary = TRUE, directions = 8, na = "outside",
               crop = FALSE)
  # Even compared against itself, an unfingerprintable source must not match.
  expect_false(isTRUE(db_cache_match(meta, meta)))
  expect_match(db_cache_match(meta, meta), "cannot be fingerprinted")
})

test_that("fingerprint = 'full' makes even a large in-memory raster cacheable", {
  skip_if_no_python()
  d <- withr::local_tempdir()
  # Over the 1e6-cell "auto" ceiling, so only "full" gives it an identity.
  r <- terra::rast(nrows = 1100, ncols = 1000, xmin = 0, xmax = 1000,
                   ymin = 0, ymax = 1100, crs = "EPSG:5070")
  set.seed(8); terra::values(r) <- rbinom(terra::ncell(r), 1, 0.3)

  analyze_patches(r, output_dir = d, metrics = FALSE, fingerprint = "auto", quiet = TRUE)
  # Under "auto" it is past the hashing ceiling, so it has no identity and the
  # cache is declined rather than guessed at.
  expect_message(analyze_patches(r, output_dir = d, metrics = FALSE,
                                 fingerprint = "auto", quiet = FALSE),
                 "cannot be fingerprinted")

  analyze_patches(r, output_dir = d, metrics = FALSE, fingerprint = "full", quiet = TRUE)
  expect_message(analyze_patches(r, output_dir = d, metrics = FALSE,
                                 fingerprint = "full", quiet = FALSE),
                 "Reusing cached result")
})
