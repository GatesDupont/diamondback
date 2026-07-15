# Temporal comparison. The design commitment here is that the overlap table is
# complete and the lineage rule is a default, not a constraint: everything
# needed to re-derive lineage under a different rule is returned.
# See DESIGN.md section 4, questions 6-8.

#' Compare patch configurations at two times
#'
#' Cross-tabulates the cells shared by patches at time 1 and time 2, identifies
#' each patch's primary predecessor and descendant, and classifies what happened
#' to it: persistence, split, merger, complex reconfiguration, appearance or
#' disappearance.
#'
#' @section What is returned, and why three tables:
#' * `$overlaps` --- one row per overlapping pair. This is the raw, complete
#'   evidence, including both directional proportions and both ranks. **No
#'   lineage rule is baked in here**, so you can re-derive lineage under any rule
#'   you prefer without re-running the expensive step.
#' * `$events` --- one row per patch (both times), with its event class, its
#'   primary predecessor or descendant, and how much of it was retained.
#' * `$lineages` --- one row per connected component of the bipartite overlap
#'   graph. A component is the honest unit for a merger: two patches becoming one
#'   is a fact about the group, not about either patch.
#'
#' @section Event classification:
#' Classes are assigned from the structure of each component, counting only links
#' that pass the thresholds:
#'
#' | class | meaning |
#' | --- | --- |
#' | `persistence` | one patch at t1, one at t2 |
#' | `split` | one patch at t1, several at t2 |
#' | `merger` | several at t1, one at t2 |
#' | `complex` | several at both --- a split and merger together |
#' | `disappearance` | at t1 with no qualifying overlap at t2 |
#' | `appearance` | at t2 with no qualifying overlap at t1 |
#'
#' `persistence` says the patches are linked, not that nothing changed: a patch
#' that lost 90% of its area but stayed connected is still `persistence`. Use
#' `prop_primary` and the area columns for magnitude.
#'
#' @section Choosing the mother patch:
#' The primary predecessor of a time-2 patch is the time-1 patch contributing the
#' most overlapping cells. Ties break deterministically --- larger time-1 patch
#' first, then lower ID --- and any broken tie is flagged in `primary_tie` so it
#' is auditable rather than silent. Primary descendants are chosen symmetrically.
#'
#' This is a default. `$overlaps` has `prop_t1_retained` and `prop_t2_inherited`
#' as well as the counts, so a rule based on proportion rather than absolute
#' overlap is a one-line re-derivation.
#'
#' @section Thresholds, and seeing what they did:
#' Defaults keep every overlap, because discarding evidence should be a decision
#' you make rather than one made for you. But classification is acutely sensitive
#' to slivers: a patch that donates 99.9% of itself to one descendant and a
#' single cell to another is called a `split` under `min_overlap_cells = 1`,
#' which is arguably not a meaningful description of what happened.
#'
#' So nothing is hidden. Thresholded links are **flagged, not deleted**:
#' `$overlaps` keeps every overlapping pair with a `passes_threshold` column.
#' And every patch is classified **twice** --- `event` uses only links that
#' passed, `event_all_overlaps` uses all of them, and
#' `threshold_changed_event` marks the patches where the two disagree.
#' `metadata$n_events_changed_by_threshold` counts them, and both thresholds are
#' stored in `metadata`.
#'
#' The point is that two analyses of the same rasters cannot disagree without the
#' cause being visible in the output. If `threshold_changed_event` is all
#' `FALSE`, your thresholds did not matter and the classification is robust. If
#' it is not, those are exactly the patches to look at.
#'
#' **Which knob suppresses a sliver.** A link passes `min_overlap_prop` if it is
#' a big enough share of **either** side, so that a small patch absorbed whole
#' into a large one is not discarded. The consequence is worth knowing: a
#' one-cell descendant is 100% of *itself*, so `min_overlap_prop` can never
#' suppress it, however large you set it. Use `min_overlap_cells` for absolute
#' slivers, and `min_overlap_prop` for links that are trivial to both sides.
#'
#' @param x1,x2 The two times. Each may be a [patch_result], or any input
#'   [label_patches()] accepts, in which case it is labelled first using
#'   `class`, `directions`, `mask` and `na`.
#' @param class,directions,mask,na Passed to [label_patches()] when `x1`/`x2`
#'   need labelling. Ignored for inputs that are already a [patch_result].
#' @param min_overlap_cells Minimum overlapping cells for a link to count.
#' @param min_overlap_prop Minimum overlap as a proportion of **either** patch
#'   for a link to count. A link survives if it is a meaningful share of one
#'   side, so a small patch absorbed into a large one is not discarded.
#' @param check_domain Warn when the two analysis domains differ, which makes
#'   appearance and disappearance ambiguous.
#' @param quiet Suppress progress reporting.
#'
#' @return A `patch_comparison` object with `$overlaps`, `$events`, `$lineages`
#'   and `$metadata`.
#'
#'   `$overlaps` columns: `id1`, `id2`, `cells1`, `cells2`, `overlap_cells`,
#'   `overlap_area_m2`, `prop_t1_retained` (share of the t1 patch now in this t2
#'   patch), `prop_t2_inherited` (share of the t2 patch that came from this t1
#'   patch), `rank_from_t1` (this t2 patch's rank among the t1 patch's
#'   descendants), `rank_from_t2` (this t1 patch's rank among the t2 patch's
#'   predecessors), `is_primary_predecessor`, `is_primary_descendant`.
#'
#'   `$events` columns: `time`, `patch_id`, `cells`, `event`, `n_links` (number
#'   of descendants for a t1 patch, of predecessors for a t2 patch),
#'   `primary_id` (mother patch for t2 rows, main descendant for t1 rows),
#'   `primary_overlap_cells`, `prop_primary` (share of this patch shared with its
#'   primary link), `primary_tie`, `lineage_id`.
#'
#' @seealso [track_patch_series()] for more than two times.
#' @export
#' @examplesIf diamondback_ready()
#' # A patch splits in two between t1 and t2.
#' m1 <- matrix(0, 5, 5); m1[2:4, 2:4] <- 1
#' m2 <- matrix(0, 5, 5); m2[2, 2:4] <- 1; m2[4, 2:4] <- 1
#' cmp <- compare_patches(m1, m2, quiet = TRUE)
#' cmp$events[, c("time", "patch_id", "event", "primary_id")]
compare_patches <- function(x1, x2,
                            class = NULL,
                            directions = 8,
                            mask = NULL,
                            na = c("outside", "background"),
                            min_overlap_cells = 1L,
                            min_overlap_prop = 0,
                            check_domain = TRUE,
                            quiet = FALSE) {
  na <- match.arg(na)
  t0 <- Sys.time()

  p1 <- db_as_patch_result(x1, "x1", class, directions, mask, na, quiet)
  p2 <- db_as_patch_result(x2, "x2", class, directions, mask, na, quiet)

  db_check_comparable(p1, p2, check_domain = check_domain, quiet = quiet)

  s1 <- db_result_state(p1, quiet = quiet)
  s2 <- db_result_state(p2, quiet = quiet)

  if (!quiet) cli::cli_alert_info("Cross-tabulating patch overlap ...")
  py <- db_py()
  ov <- py_try(
    py$overlap_counts(s1$labels, s2$labels, as.integer(s1$n), as.integer(s2$n)),
    "cross-tabulating patch overlap"
  )

  cells1 <- db_patch_cells(p1, s1)
  cells2 <- db_patch_cells(p2, s2)

  r <- s1$template
  lonlat <- db_is_lonlat(r)
  cell_area <- if (lonlat) NA_real_ else terra::xres(r) * terra::yres(r)

  overlaps <- data.frame(
    id1 = py_num(ov, "id1"),
    id2 = py_num(ov, "id2"),
    overlap_cells = py_num(ov, "cells")
  )

  n_raw <- nrow(overlaps)
  if (n_raw > 0) {
    overlaps$cells1 <- cells1[overlaps$id1]
    overlaps$cells2 <- cells2[overlaps$id2]
    overlaps$overlap_area_m2 <- overlaps$overlap_cells * cell_area
    overlaps$prop_t1_retained <- overlaps$overlap_cells / overlaps$cells1
    overlaps$prop_t2_inherited <- overlaps$overlap_cells / overlaps$cells2

    # Thresholded links are *flagged*, not deleted. The overlap table is the
    # raw evidence, and dropping rows from it would defeat that: a user could
    # not see what the threshold decided on their behalf.
    overlaps$passes_threshold <-
      overlaps$overlap_cells >= min_overlap_cells &
      (overlaps$prop_t1_retained >= min_overlap_prop |
         overlaps$prop_t2_inherited >= min_overlap_prop)
    n_dropped <- sum(!overlaps$passes_threshold)
    if (n_dropped > 0 && !quiet) {
      cli::cli_alert_info(
        "Thresholds set aside {n_dropped} of {n_raw} overlap link{?s} as incidental."
      )
    }
  } else {
    n_dropped <- 0L
    overlaps$cells1 <- numeric(0); overlaps$cells2 <- numeric(0)
    overlaps$overlap_area_m2 <- numeric(0)
    overlaps$prop_t1_retained <- numeric(0); overlaps$prop_t2_inherited <- numeric(0)
    overlaps$passes_threshold <- logical(0)
  }

  kept <- overlaps[overlaps$passes_threshold, , drop = FALSE]

  # Classify twice: once using every overlap, once using only those that passed.
  # Event classes are acutely threshold-sensitive -- one incidental cell turns
  # "persisted" into "split" -- so both readings are reported and any patch whose
  # label depends on the threshold is counted. Without this, two analyses could
  # label the same transition differently with no visible cause.
  lin_all <- db_lineages(overlaps, length(cells1), length(cells2))
  lin_kept <- db_lineages(kept, length(cells1), length(cells2))

  ranked_kept <- db_rank_overlaps(kept, cells1, cells2)
  events <- db_events(ranked_kept, cells1, cells2, lin_kept)
  events_all <- db_events(db_rank_overlaps(overlaps, cells1, cells2),
                          cells1, cells2, lin_all)

  events$event_all_overlaps <- events_all$event
  events$n_links_all_overlaps <- events_all$n_links
  events$threshold_changed_event <- events$event != events$event_all_overlaps
  n_changed <- sum(events$threshold_changed_event)

  # Carry the ranks back onto the full table; non-passing links have no rank
  # because lineage is not built from them.
  overlaps <- db_merge_ranks(overlaps, ranked_kept)

  ord <- order(overlaps$id1, -overlaps$overlap_cells)
  overlaps <- overlaps[ord, , drop = FALSE]
  rownames(overlaps) <- NULL
  lineage <- lin_kept

  meta <- list(
    created = as.character(Sys.time()),
    package_version = as.character(utils::packageVersion("diamondback")),
    elapsed_secs = db_elapsed_num(t0),
    n_patches_t1 = length(cells1),
    n_patches_t2 = length(cells2),
    n_links = nrow(kept),
    n_links_all = n_raw,
    n_links_dropped = n_dropped,
    # The knobs that produced these labels, stored with them. A comparison
    # should never be reproducible only by remembering what you typed.
    min_overlap_cells = min_overlap_cells,
    min_overlap_prop = min_overlap_prop,
    n_events_changed_by_threshold = n_changed,
    geometry = db_geometry(r),
    t1 = .subset2(p1, "metadata"),
    t2 = .subset2(p2, "metadata")
  )

  res <- structure(
    list(overlaps = overlaps, events = events, lineages = lineage$table,
         metadata = meta),
    class = "patch_comparison"
  )
  if (!quiet) {
    tb <- table(events$event)
    cli::cli_alert_success("Compared in {db_elapsed(t0)}: {paste(names(tb), as.integer(tb), sep = ' = ', collapse = ', ')}.")
  }
  res
}

# ---------------------------------------------------------------------------
# internals
# ---------------------------------------------------------------------------

db_as_patch_result <- function(x, arg, class, directions, mask, na, quiet) {
  if (inherits(x, "patch_result")) return(x)
  if (!quiet) cli::cli_alert_info("Labelling {.arg {arg}} ...")
  label_patches(x, class = class, directions = directions, mask = mask,
                na = na, crop = FALSE, quiet = quiet)
}

#' Refuse to compare rasters that do not line up
#'
#' Silent resampling is the failure mode this exists to prevent: two rasters
#' that differ by half a cell will still produce a full, plausible, wrong
#' overlap table.
#' @noRd
db_check_comparable <- function(p1, p2, check_domain = TRUE, quiet = FALSE) {
  g1 <- .subset2(p1, "metadata")$geometry
  g2 <- .subset2(p2, "metadata")$geometry
  d <- db_geometry_diff(g1, g2)
  if (length(d)) {
    cli::cli_abort(c(
      "The two rasters are not aligned, so their cells cannot be compared.",
      "x" = "Mismatched: {.val {d}}.",
      "i" = "diamondback will not resample silently: a half-cell offset would still \\
             produce a complete and entirely wrong overlap table.",
      "*" = "Align them first, e.g. {.code x2 <- terra::resample(x2, x1, method = \"near\")}.",
      "*" = "Or crop both to a common extent with {.fn terra::crop}."
    ), call = NULL)
  }

  m1 <- .subset2(p1, "metadata")
  m2 <- .subset2(p2, "metadata")
  if (!identical(m1$directions, m2$directions) && !is.na(m1$directions) && !is.na(m2$directions)) {
    cli::cli_warn(c(
      "The two results used different connectivity ({m1$directions} vs {m2$directions}).",
      "i" = "Patch identity is not comparable across connectivity rules."
    ))
  }

  if (isTRUE(check_domain)) {
    n1 <- m1$cells$outside + m1$cells$missing
    n2 <- m2$cells$outside + m2$cells$missing
    if (!identical(n1, n2)) {
      cli::cli_warn(c(
        "The two analysis domains differ: {format(n1, big.mark = ',')} vs \\
         {format(n2, big.mark = ',')} cells excluded.",
        "i" = "A patch can then appear or disappear purely because the domain changed.",
        "i" = "Use the same {.arg mask} for both times, or set {.code check_domain = FALSE} \\
               if the difference is intended."
      ))
    }
  }
  invisible(TRUE)
}

db_patch_cells <- function(p, st) {
  met <- .subset2(p, "metrics")
  if (!is.null(met) && "cells" %in% names(met) && nrow(met) == st$n && !anyNA(met$cells)) {
    return(as.numeric(met$cells))
  }
  py_num(py_try(db_py()$patch_stats(st$labels, as.integer(st$n)),
                "counting patch cells"), "count")[-1]
}

#' Put ranks from the passing links back onto the full overlap table
#'
#' Non-passing links get NA ranks rather than a rank of their own: lineage is
#' built only from links that passed, so a rank for a link that took no part in
#' it would be a fiction.
#' @noRd
db_merge_ranks <- function(all_ov, ranked) {
  n <- nrow(all_ov)
  int_cols <- c("rank_from_t1", "rank_from_t2")
  lgl_cols <- c("is_primary_descendant", "is_primary_predecessor",
                "tie_from_t1", "tie_from_t2")
  # rep() rather than a scalar: a zero-row overlap table (nothing overlapped at
  # all) must still come back with the full set of columns.
  for (cn in int_cols) all_ov[[cn]] <- rep(NA_integer_, n)
  for (cn in lgl_cols) all_ov[[cn]] <- rep(NA, n)
  if (!nrow(ranked) || !n) return(all_ov)

  key_all <- paste(all_ov$id1, all_ov$id2, sep = "_")
  key_rk <- paste(ranked$id1, ranked$id2, sep = "_")
  i <- match(key_rk, key_all)
  for (cn in c(int_cols, lgl_cols)) all_ov[[cn]][i] <- ranked[[cn]]
  all_ov
}

#' Rank descendants and predecessors, and mark the primary of each
#'
#' Ties are broken by patch size then ID so that results are reproducible.
#' `primary_tie` records that a tie-break was needed.
#' @noRd
db_rank_overlaps <- function(ov, cells1, cells2) {
  if (!nrow(ov)) {
    ov$rank_from_t1 <- integer(0); ov$rank_from_t2 <- integer(0)
    ov$is_primary_descendant <- logical(0); ov$is_primary_predecessor <- logical(0)
    ov$tie_from_t1 <- logical(0); ov$tie_from_t2 <- logical(0)
    return(ov)
  }

  # Rank this t2 patch among the descendants of its t1 patch.
  ov$rank_from_t1 <- db_rank_within(ov$id1, ov$overlap_cells, cells2[ov$id2], ov$id2)
  # Rank this t1 patch among the predecessors of its t2 patch.
  ov$rank_from_t2 <- db_rank_within(ov$id2, ov$overlap_cells, cells1[ov$id1], ov$id1)

  ov$is_primary_descendant <- ov$rank_from_t1 == 1L
  ov$is_primary_predecessor <- ov$rank_from_t2 == 1L

  ov$tie_from_t1 <- db_tie_flag(ov$id1, ov$overlap_cells)
  ov$tie_from_t2 <- db_tie_flag(ov$id2, ov$overlap_cells)
  ov
}

# Rank rows within groups by overlap desc, then size desc, then id asc.
db_rank_within <- function(group, overlap, size, id) {
  o <- order(group, -overlap, -size, id)
  r <- integer(length(group))
  g <- group[o]
  # Position within each run of equal group values.
  runs <- rle(g)
  r[o] <- sequence(runs$lengths)
  r
}

# TRUE when a row's top overlap within its group is tied with another row.
db_tie_flag <- function(group, overlap) {
  mx <- stats::ave(overlap, group, FUN = max)
  n_at_max <- stats::ave(as.numeric(overlap == mx), group, FUN = sum)
  (overlap == mx) & (n_at_max > 1)
}

#' Connected components of the bipartite overlap graph
#'
#' Union-find over t1 and t2 patch nodes. Patches with no qualifying link get a
#' component of their own, which is what makes appearance and disappearance fall
#' out of the same machinery as everything else.
#' @noRd
db_lineages <- function(ov, n1, n2) {
  total <- n1 + n2
  if (total == 0L) {
    return(list(comp1 = integer(0), comp2 = integer(0),
                table = data.frame(lineage_id = integer(), n_t1 = integer(),
                                   n_t2 = integer(), event = character(),
                                   cells_t1 = numeric(), cells_t2 = numeric())))
  }

  parent <- seq_len(total)
  find <- function(i) {
    root <- as.integer(i)
    while (parent[root] != root) root <- parent[root]
    while (parent[i] != root) {
      nxt <- parent[i]
      parent[i] <<- root
      i <- nxt
    }
    root
  }
  union <- function(a, b) {
    ra <- find(a); rb <- find(b)
    # Patch IDs arrive as doubles from NumPy; keep `parent` integer throughout,
    # or path compression silently promotes it and vapply's type check fails.
    if (ra != rb) parent[rb] <<- as.integer(ra)
    invisible(NULL)
  }

  if (nrow(ov)) {
    for (k in seq_len(nrow(ov))) {
      union(as.integer(ov$id1[k]), as.integer(n1 + ov$id2[k]))
    }
  }

  roots <- vapply(seq_len(total), find, integer(1))
  comp <- match(roots, sort(unique(roots)))

  comp1 <- if (n1) comp[seq_len(n1)] else integer(0)
  comp2 <- if (n2) comp[(n1 + 1L):total] else integer(0)

  n_comp <- max(comp)
  n_t1 <- tabulate(comp1, nbins = n_comp)
  n_t2 <- tabulate(comp2, nbins = n_comp)

  event <- rep(NA_character_, n_comp)
  event[n_t1 == 1 & n_t2 == 1] <- "persistence"
  event[n_t1 == 1 & n_t2 > 1] <- "split"
  event[n_t1 > 1 & n_t2 == 1] <- "merger"
  event[n_t1 > 1 & n_t2 > 1] <- "complex"
  event[n_t1 >= 1 & n_t2 == 0] <- "disappearance"
  event[n_t1 == 0 & n_t2 >= 1] <- "appearance"

  list(
    comp1 = comp1, comp2 = comp2,
    table = data.frame(
      lineage_id = seq_len(n_comp),
      n_t1 = n_t1,
      n_t2 = n_t2,
      event = event
    )
  )
}

#' One row per patch, at both times
#' @noRd
db_events <- function(ov, cells1, cells2, lineage) {
  n1 <- length(cells1); n2 <- length(cells2)
  ev_by_comp <- lineage$table$event

  mk <- function(time, cells, comp, is_t1) {
    n <- length(cells)
    if (n == 0L) {
      return(data.frame(time = character(), patch_id = integer(), cells = numeric(),
                        event = character(), n_links = integer(), primary_id = integer(),
                        primary_overlap_cells = numeric(), prop_primary = numeric(),
                        primary_tie = logical(), lineage_id = integer()))
    }
    prim_id <- rep(NA_integer_, n)
    prim_ov <- rep(NA_real_, n)
    prim_tie <- rep(FALSE, n)
    n_links <- integer(n)

    if (nrow(ov)) {
      self <- if (is_t1) ov$id1 else ov$id2
      other <- if (is_t1) ov$id2 else ov$id1
      # For a t1 patch the "primary" is its main descendant; for a t2 patch it is
      # its mother patch. Both are rank 1 in the corresponding direction.
      rank <- if (is_t1) ov$rank_from_t1 else ov$rank_from_t2
      tie <- if (is_t1) ov$tie_from_t1 else ov$tie_from_t2

      n_links <- tabulate(self, nbins = n)
      top <- rank == 1L
      prim_id[self[top]] <- as.integer(other[top])
      prim_ov[self[top]] <- ov$overlap_cells[top]
      prim_tie[self[top]] <- tie[top]
    }

    data.frame(
      time = time,
      patch_id = seq_len(n),
      cells = cells,
      event = ev_by_comp[comp],
      n_links = n_links,
      primary_id = prim_id,
      primary_overlap_cells = prim_ov,
      prop_primary = prim_ov / cells,
      primary_tie = prim_tie,
      lineage_id = comp
    )
  }

  rbind(
    mk("t1", cells1, lineage$comp1, TRUE),
    mk("t2", cells2, lineage$comp2, FALSE)
  )
}

#' @export
print.patch_comparison <- function(x, ...) {
  m <- .subset2(x, "metadata")
  cli::cli_h1("<patch_comparison>")
  cli::cli_dl(c(
    "time 1" = "{format(m$n_patches_t1, big.mark = ',')} patches",
    "time 2" = "{format(m$n_patches_t2, big.mark = ',')} patches",
    "overlap links" = "{format(m$n_links, big.mark = ',')}\\
                       {if (m$n_links_dropped > 0) paste0(' (', m$n_links_dropped, ' dropped by thresholds)') else ''}"
  ))

  ev <- .subset2(x, "events")
  if (nrow(ev)) {
    cli::cli_h2("Events")
    tb <- table(factor(ev$event, levels = c("persistence", "split", "merger",
                                            "complex", "appearance", "disappearance")))
    tb <- tb[tb > 0]
    cli::cli_dl(stats::setNames(as.character(as.integer(tb)), names(tb)))
    n_tie <- sum(ev$primary_tie, na.rm = TRUE)
    if (n_tie > 0) {
      cli::cli_alert_warning("{n_tie} patch{?es} had a tied primary link, resolved by size then ID.")
    }
    n_chg <- m$n_events_changed_by_threshold %||% 0L
    if (n_chg > 0) {
      cli::cli_alert_warning(
        "{n_chg} patch{?es} would be classified differently without the overlap \\
         thresholds (see {.code $events$event_all_overlaps})."
      )
    }
  }
  cli::cli_text("")
  cli::cli_text("{.strong $overlaps} {nrow(.subset2(x, 'overlaps'))} x {ncol(.subset2(x, 'overlaps'))} \\
                 {.emph |} {.strong $events} {nrow(ev)} x {ncol(ev)} \\
                 {.emph |} {.strong $lineages} {nrow(.subset2(x, 'lineages'))} x {ncol(.subset2(x, 'lineages'))}")
  invisible(x)
}

#' @export
summary.patch_comparison <- function(object, ...) {
  print(object)
  ev <- .subset2(object, "events")
  if (!nrow(ev)) return(invisible(object))

  cli::cli_h2("Area accounting")
  c1 <- sum(ev$cells[ev$time == "t1"])
  c2 <- sum(ev$cells[ev$time == "t2"])
  cli::cli_dl(c(
    "patch cells at t1" = format(c1, big.mark = ","),
    "patch cells at t2" = format(c2, big.mark = ","),
    "net change" = sprintf("%+s (%+.1f%%)", format(c2 - c1, big.mark = ","),
                           100 * (c2 - c1) / max(c1, 1))
  ))
  t2 <- ev[ev$time == "t2", ]
  if (any(!is.na(t2$prop_primary))) {
    cli::cli_dl(c(
      "median share of a t2 patch inherited from its mother" =
        sprintf("%.2f", stats::median(t2$prop_primary, na.rm = TRUE))
    ))
  }
  invisible(object)
}
