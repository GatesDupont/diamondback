# The tiled distance transform must be *identical* to the full-array one, not
# close to it. Everything here forces absurdly small tiles so that ordinary test
# rasters cross many seams, then demands the same answer as an untiled run.
#
# `target_bytes` is the lever: it sets the strip working set, so a tiny value
# forces one-row strips out of a 60-row raster.

core_counts_at <- function(st, depth, edge = "all", sampling = c(1, 1),
                           target_bytes = 1e9) {
  py <- diamondback:::db_py()
  out <- py$core_counts(st$code, st$labels, as.integer(st$n), 0L, depth,
                        reticulate::tuple(sampling[1], sampling[2]), edge,
                        want_mask = FALSE, target_bytes = as.integer(target_bytes))
  py_num(out, "core_count")
}

state_of <- function(m, res = 1, yres = res) {
  diamondback:::db_result_state(
    label_patches(rast_from(m, res = res, yres = yres), quiet = TRUE), quiet = TRUE)
}

test_that("tiling reproduces the full-array transform exactly", {
  skip_if_no_python()
  set.seed(11)
  for (trial in 1:4) {
    m <- matrix(rbinom(60 * 50, 1, runif(1, 0.3, 0.8)), 60, 50)
    st <- state_of(m)
    for (depth in c(0.5, 1, 2, 3.5)) {
      for (edge in c("all", "background")) {
        untiled <- core_counts_at(st, depth, edge, target_bytes = 1e9)
        tiled   <- core_counts_at(st, depth, edge, target_bytes = 2000)
        expect_identical(untiled, tiled,
                         info = sprintf("trial %d, depth %s, edge %s", trial, depth, edge))
      }
    }
  }
})

test_that("tiling is exact with anisotropic cells", {
  skip_if_no_python()
  # The halo is derived from the row spacing, so a tall thin cell is where an
  # off-by-one in that derivation would show up.
  set.seed(12)
  m <- matrix(rbinom(60 * 40, 1, 0.6), 60, 40)
  st <- state_of(m, res = 1, yres = 10)
  for (depth in c(1, 5, 12, 25)) {
    for (samp in list(c(10, 1), c(1, 10), c(3, 7))) {
      expect_identical(
        core_counts_at(st, depth, "all", samp, target_bytes = 1e9),
        core_counts_at(st, depth, "all", samp, target_bytes = 1500),
        info = sprintf("depth %s, sampling %s", depth, paste(samp, collapse = "x"))
      )
    }
  }
})

test_that("tiling is exact when the halo is larger than the tile", {
  skip_if_no_python()
  # A deep edge depth against a tiny tile budget: the halo then dominates the
  # window, which is the degenerate case for the strip arithmetic.
  set.seed(13)
  m <- matrix(rbinom(50 * 50, 1, 0.75), 50, 50)
  st <- state_of(m)
  for (depth in c(8, 15, 30)) {
    expect_identical(
      core_counts_at(st, depth, "all", target_bytes = 1e9),
      core_counts_at(st, depth, "all", target_bytes = 500)
    )
  }
})

test_that("the grid border is still an edge source under tiling", {
  skip_if_no_python()
  # Only true borders may be padded. If interior strip seams were padded too,
  # cells beside a seam would look edge-affected and core area would collapse.
  m <- matrix(1, 40, 40)          # solid habitat: every edge is the grid border
  st <- state_of(m)

  untiled <- core_counts_at(st, 3, "all", target_bytes = 1e9)
  tiled   <- core_counts_at(st, 3, "all", target_bytes = 800)
  expect_identical(untiled, tiled)

  # A 40x40 block eroded by >3 leaves the 34x34 interior. If seams were padded,
  # the tiled count would be far lower.
  expect_equal(unname(tiled[2]), 34 * 34)
})

test_that("a patch spanning many strips keeps its core intact", {
  skip_if_no_python()
  # One tall column of habitat crossing every seam.
  m <- matrix(0, 80, 11); m[, 2:10] <- 1
  st <- state_of(m)
  untiled <- core_counts_at(st, 2, "all", target_bytes = 1e9)
  tiled   <- core_counts_at(st, 2, "all", target_bytes = 400)
  expect_identical(untiled, tiled)
  expect_gt(unname(tiled[2]), 0)
})

test_that("edge = 'background' tiles exactly with a mask present", {
  skip_if_no_python()
  set.seed(14)
  v <- matrix(rbinom(60 * 40, 1, 0.7), 60, 40)
  mk <- matrix(1, 60, 40); mk[1:10, ] <- 0; mk[, 35:40] <- 0
  res <- label_patches(rast_from(v), mask = rast_from(mk), crop = FALSE, quiet = TRUE)
  st <- diamondback:::db_result_state(res, quiet = TRUE)
  for (edge in c("all", "background")) {
    expect_identical(
      core_counts_at(st, 2, edge, target_bytes = 1e9),
      core_counts_at(st, 2, edge, target_bytes = 900)
    )
  }
})

test_that("the core raster matches the counts under tiling", {
  skip_if_no_python()
  m <- matrix(0, 40, 40); m[5:35, 5:35] <- 1
  res <- patch_core_area(label_patches(rast_from(m), quiet = TRUE), depth = 2,
                         core_raster = TRUE, quiet = TRUE)
  v <- terra::values(res$core_raster, mat = FALSE)
  expect_equal(sum(v == 1), res$metrics$core_cells[1])
})

test_that("the memory estimate reflects the tiled transform", {
  # The full-array float64 distance grid is gone, so core is no longer the most
  # expensive stage by a factor of two.
  expect_equal(db_estimate_memory(1e6, "core"), 1e6 * 6)
  expect_lt(db_estimate_memory(1e6, "core"), 1e6 * 13)   # the old figure
  expect_equal(db_estimate_memory(1e6, "core"), db_estimate_memory(1e6, "label"))
})

test_that("core area at 300 m on a Central-Hardwoods-sized raster is not refused", {
  # The regression that started this: 724M cells needed 8.8 GB and could not run
  # on an 8 GB machine. It must now pass the guard on that machine.
  need <- db_estimate_memory(724461424, "core")
  expect_lt(need, 8 * 1024^3 * 0.9)
  expect_silent(db_check_memory(724461424, "core", memory_limit = 8 * 1024^3, quiet = TRUE))
})
