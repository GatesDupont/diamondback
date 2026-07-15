test_that("an unchanged patch is persistence with full retention", {
  skip_if_no_python()
  m <- matrix(0, 5, 5); m[2:4, 2:4] <- 1
  cmp <- compare_patches(m, m, quiet = TRUE)

  expect_equal(nrow(cmp$overlaps), 1L)
  expect_equal(cmp$overlaps$overlap_cells, 9)
  expect_equal(cmp$overlaps$prop_t1_retained, 1)
  expect_equal(cmp$overlaps$prop_t2_inherited, 1)
  expect_true(all(cmp$events$event == "persistence"))
  expect_equal(nrow(cmp$lineages), 1L)
})

test_that("a shrinking but connected patch is still persistence", {
  skip_if_no_python()
  m1 <- matrix(0, 6, 6); m1[2:5, 2:5] <- 1
  m2 <- matrix(0, 6, 6); m2[2:3, 2:3] <- 1
  cmp <- compare_patches(m1, m2, quiet = TRUE)
  expect_true(all(cmp$events$event == "persistence"))
  expect_equal(cmp$overlaps$prop_t1_retained, 4 / 16)
  expect_equal(cmp$overlaps$prop_t2_inherited, 1)
})

test_that("a split is classified and both descendants trace to the mother", {
  skip_if_no_python()
  m1 <- matrix(0, 5, 5); m1[2:4, 2:4] <- 1
  m2 <- matrix(0, 5, 5); m2[2, 2:4] <- 1; m2[4, 2:4] <- 1
  cmp <- compare_patches(m1, m2, quiet = TRUE)

  expect_true(all(cmp$events$event == "split"))
  expect_equal(nrow(cmp$overlaps), 2L)

  t1 <- cmp$events[cmp$events$time == "t1", ]
  expect_equal(t1$n_links, 2L)          # two descendants

  t2 <- cmp$events[cmp$events$time == "t2", ]
  expect_equal(nrow(t2), 2L)
  expect_true(all(t2$primary_id == 1))  # both descend from t1 patch 1
  expect_true(all(t2$prop_primary == 1))
})

test_that("a merger is classified and both predecessors share a descendant", {
  skip_if_no_python()
  m1 <- matrix(0, 5, 5); m1[2, 2:4] <- 1; m1[4, 2:4] <- 1
  m2 <- matrix(0, 5, 5); m2[2:4, 2:4] <- 1
  cmp <- compare_patches(m1, m2, quiet = TRUE)

  expect_true(all(cmp$events$event == "merger"))
  t2 <- cmp$events[cmp$events$time == "t2", ]
  expect_equal(t2$n_links, 2L)          # two predecessors
  expect_equal(nrow(cmp$lineages), 1L)
  expect_equal(cmp$lineages$n_t1, 2L)
  expect_equal(cmp$lineages$n_t2, 1L)
})

test_that("a simultaneous split and merger is classified as complex", {
  skip_if_no_python()
  # Two t1 bars; at t2 they are joined by a bridge and also split off a piece,
  # leaving one component with 2 t1 patches and 2 t2 patches.
  m1 <- matrix(0, 7, 7)
  m1[2, 2:6] <- 1        # t1 patch A (top bar)
  m1[6, 2:6] <- 1        # t1 patch B (bottom bar)

  m2 <- matrix(0, 7, 7)
  m2[2, 2:3] <- 1; m2[3:6, 3] <- 1; m2[6, 2:3] <- 1  # left: joins A and B
  m2[2, 5:6] <- 1; m2[3:6, 5] <- 1; m2[6, 5:6] <- 1  # right: also joins A and B

  cmp <- compare_patches(m1, m2, quiet = TRUE)
  expect_equal(nrow(cmp$lineages), 1L)
  expect_equal(cmp$lineages$n_t1, 2L)
  expect_equal(cmp$lineages$n_t2, 2L)
  expect_equal(cmp$lineages$event, "complex")
  expect_true(all(cmp$events$event == "complex"))
})

test_that("disappearance and appearance are detected when nothing overlaps", {
  skip_if_no_python()
  m1 <- matrix(0, 5, 5); m1[1:2, 1:2] <- 1
  m2 <- matrix(0, 5, 5); m2[4:5, 4:5] <- 1
  cmp <- compare_patches(m1, m2, quiet = TRUE)

  expect_equal(nrow(cmp$overlaps), 0L)
  expect_equal(cmp$events$event[cmp$events$time == "t1"], "disappearance")
  expect_equal(cmp$events$event[cmp$events$time == "t2"], "appearance")
  expect_true(is.na(cmp$events$primary_id[1]))
  expect_equal(cmp$events$n_links, c(0L, 0L))
})

test_that("complete disappearance leaves only t1 rows", {
  skip_if_no_python()
  m1 <- matrix(0, 4, 4); m1[2:3, 2:3] <- 1
  m2 <- matrix(0, 4, 4)
  expect_warning(cmp <- compare_patches(m1, m2, quiet = TRUE), "No foreground")
  expect_equal(nrow(cmp$events), 1L)
  expect_equal(cmp$events$event, "disappearance")
})

test_that("appearance from nothing leaves only t2 rows", {
  skip_if_no_python()
  m1 <- matrix(0, 4, 4)
  m2 <- matrix(0, 4, 4); m2[2:3, 2:3] <- 1
  expect_warning(cmp <- compare_patches(m1, m2, quiet = TRUE), "No foreground")
  expect_equal(nrow(cmp$events), 1L)
  expect_equal(cmp$events$event, "appearance")
})

# --- mother patch selection -------------------------------------------------

test_that("the mother patch is the largest contributor, not the largest patch", {
  skip_if_no_python()
  # t1: a big patch (rows 1-2) and a small one (row 4).
  # t2: one patch overlapping the small one by 3 cells and the big one by 1.
  m1 <- matrix(0, 5, 5)
  m1[1:2, 1:4] <- 1     # big: 8 cells
  m1[4, 1:3] <- 1       # small: 3 cells
  m2 <- matrix(0, 5, 5)
  m2[2, 1] <- 1         # 1 cell overlapping the big patch
  m2[3, 1] <- 1         # bridge
  m2[4, 1:3] <- 1       # 3 cells overlapping the small patch

  cmp <- compare_patches(m1, m2, directions = 4, quiet = TRUE)
  t2 <- cmp$events[cmp$events$time == "t2", ]
  ov <- cmp$overlaps

  big_id <- 1; small_id <- 2
  expect_equal(ov$overlap_cells[ov$id1 == big_id], 1)
  expect_equal(ov$overlap_cells[ov$id1 == small_id], 3)
  # The mother is the small patch: it contributed more cells.
  expect_equal(t2$primary_id, small_id)
})

test_that("ties in overlap are broken deterministically and flagged", {
  skip_if_no_python()
  # A perfectly symmetric split: both descendants get 3 cells.
  m1 <- matrix(0, 5, 5); m1[2:4, 2:4] <- 1
  m2 <- matrix(0, 5, 5); m2[2, 2:4] <- 1; m2[4, 2:4] <- 1
  cmp <- compare_patches(m1, m2, quiet = TRUE)

  t1 <- cmp$events[cmp$events$time == "t1", ]
  expect_true(t1$primary_tie)             # the choice of main descendant was a tie
  expect_false(is.na(t1$primary_id))      # but it still made a definite choice

  # Deterministic: the same inputs give the same answer every time.
  again <- compare_patches(m1, m2, quiet = TRUE)
  expect_equal(cmp$events$primary_id, again$events$primary_id)
})

test_that("overlaps carry both directional proportions and both ranks", {
  skip_if_no_python()
  m1 <- matrix(0, 5, 5); m1[2:4, 2:4] <- 1
  m2 <- matrix(0, 5, 5); m2[2, 2:4] <- 1; m2[4, 2:4] <- 1
  ov <- compare_patches(m1, m2, quiet = TRUE)$overlaps

  expect_true(all(c("prop_t1_retained", "prop_t2_inherited",
                    "rank_from_t1", "rank_from_t2",
                    "is_primary_predecessor", "is_primary_descendant") %in% names(ov)))
  expect_equal(sort(ov$rank_from_t1), 1:2)      # ranked among t1's descendants
  expect_true(all(ov$rank_from_t2 == 1))        # each t2 patch has one predecessor
  expect_equal(sum(ov$is_primary_descendant), 1)
})

# --- thresholds -------------------------------------------------------------

test_that("min_overlap_cells suppresses trivial links and changes the event", {
  skip_if_no_python()
  # t1: one 9-cell patch at rows 2-4, cols 2-4.
  # t2: patch A (rows 2-3, cols 2-3) overlapping it by 4 cells, and patch B,
  # a single cell at (4,4) inside t1 but only diagonally touching A, so under
  # 4-connectivity it is a separate patch overlapping t1 by exactly 1 cell.
  m1 <- matrix(0, 6, 6); m1[2:4, 2:4] <- 1
  m2 <- matrix(0, 6, 6); m2[2:3, 2:3] <- 1; m2[4, 4] <- 1

  loose <- compare_patches(m1, m2, directions = 4, min_overlap_cells = 1, quiet = TRUE)
  strict <- compare_patches(m1, m2, directions = 4, min_overlap_cells = 2, quiet = TRUE)

  # Both tables keep every overlapping pair: thresholds flag, they do not delete.
  expect_equal(nrow(loose$overlaps), 2L)
  expect_equal(nrow(strict$overlaps), 2L)
  expect_equal(sort(loose$overlaps$overlap_cells), c(1, 4))
  expect_true(all(loose$overlaps$passes_threshold))
  expect_equal(sum(strict$overlaps$passes_threshold), 1L)
  expect_equal(strict$metadata$n_links_dropped, 1L)
  expect_equal(strict$metadata$n_links_all, 2L)
  expect_equal(strict$metadata$n_links, 1L)

  # And the effect of the threshold is visible per patch, not just in a count.
  expect_true(all(strict$events$threshold_changed_event))
  expect_true(all(strict$events$event_all_overlaps == "split"))
  expect_gt(strict$metadata$n_events_changed_by_threshold, 0L)

  # This is the sensitivity the threshold exists for: one incidental cell is
  # the difference between "the patch split" and "the patch persisted".
  expect_true("split" %in% loose$events$event)
  expect_false("split" %in% strict$events$event)
  expect_true("persistence" %in% strict$events$event)
  expect_true("appearance" %in% strict$events$event)
})

test_that("min_overlap_prop keeps a link that matters to either side", {
  skip_if_no_python()
  m1 <- matrix(0, 6, 6); m1[2:5, 2:5] <- 1     # 16 cells
  m2 <- matrix(0, 6, 6); m2[5, 5] <- 1         # 1 cell, fully inside t1

  # The link is 1/16 of t1 but 1/1 of t2, so a 50% threshold keeps it.
  cmp <- compare_patches(m1, m2, min_overlap_prop = 0.5, quiet = TRUE)
  expect_equal(nrow(cmp$overlaps), 1L)
  expect_equal(cmp$overlaps$prop_t2_inherited, 1)
})

test_that("thresholds report how many links they set aside", {
  skip_if_no_python()
  m1 <- matrix(0, 6, 6); m1[2:5, 2:5] <- 1
  m2 <- matrix(0, 6, 6); m2[5, 5] <- 1
  cmp <- compare_patches(m1, m2, min_overlap_prop = 1.5, quiet = TRUE)
  expect_equal(cmp$metadata$n_links_dropped, 1L)
  expect_equal(cmp$metadata$n_links, 0L)
  # The evidence is still there; only the flag says it was not used.
  expect_equal(nrow(cmp$overlaps), 1L)
  expect_false(cmp$overlaps$passes_threshold)
  # A link that took no part in lineage gets no rank rather than a fictional one.
  expect_true(is.na(cmp$overlaps$rank_from_t1))
})

# --- validation of inputs ---------------------------------------------------

test_that("misaligned rasters error instead of being resampled", {
  skip_if_no_python()
  r1 <- rast_from(matrix(1, 4, 4))
  r2 <- terra::rast(nrows = 4, ncols = 4, xmin = 5, xmax = 45, ymin = 0, ymax = 40,
                    crs = "EPSG:5070")
  terra::values(r2) <- 1
  expect_error(compare_patches(r1, r2, quiet = TRUE), "not aligned")
})

test_that("different resolutions error", {
  skip_if_no_python()
  r1 <- rast_from(matrix(1, 4, 4), res = 10)
  r2 <- rast_from(matrix(1, 4, 4), res = 20)
  expect_error(compare_patches(r1, r2, quiet = TRUE), "not aligned")
})

test_that("different dimensions error", {
  skip_if_no_python()
  expect_error(compare_patches(matrix(1, 4, 4), matrix(1, 5, 5), quiet = TRUE),
               "not aligned")
})

test_that("differing analysis domains warn", {
  skip_if_no_python()
  v <- rast_from(matrix(1, 4, 4))
  mk1 <- rast_from(matrix(1, 4, 4))
  mk2 <- rast_from(matrix(c(0, rep(1, 15)), 4, 4))

  p1 <- label_patches(v, mask = mk1, crop = FALSE, quiet = TRUE)
  p2 <- label_patches(v, mask = mk2, crop = FALSE, quiet = TRUE)
  expect_warning(compare_patches(p1, p2, quiet = TRUE), "domains differ")
})

test_that("compare_patches accepts patch_result objects directly", {
  skip_if_no_python()
  m1 <- matrix(0, 5, 5); m1[2:4, 2:4] <- 1
  m2 <- matrix(0, 5, 5); m2[2:3, 2:3] <- 1
  p1 <- analyze_patches(m1, quiet = TRUE)
  p2 <- analyze_patches(m2, quiet = TRUE)
  cmp <- compare_patches(p1, p2, quiet = TRUE)
  expect_s3_class(cmp, "patch_comparison")
  expect_equal(nrow(cmp$overlaps), 1L)
})

test_that("overlap areas use cell size", {
  skip_if_no_python()
  m <- matrix(0, 5, 5); m[2:4, 2:4] <- 1
  cmp <- compare_patches(rast_from(m, res = 10), rast_from(m, res = 10), quiet = TRUE)
  expect_equal(cmp$overlaps$overlap_area_m2, 900)   # 9 cells x 100 m2
})
