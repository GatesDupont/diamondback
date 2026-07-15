# Connected-component labelling. This is the expensive step and the one worth
# caching; see DESIGN.md sections 3 and 6.

#' Label contiguous patches in a categorical or binary raster
#'
#' Connected-component labelling via `scipy.ndimage.label()`, with explicit
#' handling of missing values and study-area masks. This is the workhorse of the
#' package: on rasters of hundreds of millions of cells it is typically far
#' faster than [terra::patches()], and it produces the same partition.
#'
#' @section Cell states:
#' diamondback keeps four cell states apart, and never lets `NA` silently become
#' foreground or background:
#'
#' * **patch** --- a foreground cell, labelled `1..N`;
#' * **background** --- a valid, in-domain cell that is not the target class, `0`;
#' * **outside the domain** --- excluded by `mask`, `NA` in the returned raster;
#' * **missing** --- `NA` in the source but inside the mask, also `NA` in the
#'   returned raster.
#'
#' The returned raster cannot distinguish the last two, because a raster has one
#' no-data value. Exact counts for all four are kept in `result$metadata$cells`,
#' and [patch_domain()] reconstructs the distinction as a categorical raster.
#'
#' `NA` is never foreground. It is not background either unless you say so with
#' `na = "background"`, which is the right choice when a raster uses `NA` to mean
#' genuine absence rather than lack of knowledge.
#'
#' @section Patch IDs are labels, not identities:
#' IDs are assigned in raster scan order (top-left to bottom-right), so they are
#' **deterministic** for a given raster, class, mask and connectivity: the same
#' inputs always produce the same numbering.
#'
#' They are **not** a persistent identity. Patch 47 in 1985 and patch 47 in 2024
#' are unrelated. Worse, IDs are not even stable across changes you might think
#' are cosmetic: cropping the raster, changing the mask, switching connectivity,
#' or a single cell flipping near the top-left will renumber everything after it.
#'
#' Never join two runs on `patch_id`. To follow patches through time, use
#' [compare_patches()] or [track_patch_series()], which establish correspondence
#' from actual cell overlap. `lineage_id` in a comparison identifies a group of
#' related patches *within that comparison*, and is likewise not a durable key.
#'
#' @section Determinism and threading:
#' Every kernel is single-threaded and deterministic. `scipy.ndimage.label()`,
#' the distance transform, and the NumPy reductions used here are serial C code,
#' and none of them touch BLAS, so `OMP_NUM_THREADS` and friends have no effect
#' on the result or the runtime. Re-running an analysis on the same inputs gives
#' bit-for-bit identical output, on any machine with a compatible NumPy and
#' SciPy. Version 1 uses one process and one core.
#'
#' @param x A `SpatRaster`, a raster filename, or a matrix.
#' @param class Which values are foreground. **One element of `class` is one
#'   class**:
#'
#'   * `NULL` (default) --- treat `x` as binary: any non-zero, non-`NA` value is
#'     foreground.
#'   * a numeric vector --- **one** class made of those values, e.g.
#'     `class = c(41, 42, 43)` labels deciduous, evergreen and mixed forest as a
#'     single forest class, connected across the boundaries between them.
#'   * a list --- several classes, labelled separately and never joined to each
#'     other, e.g. `class = list(forest = c(41, 42, 43), wetland = c(90, 95))`.
#'     Names become the `class` column in `$metrics`; unnamed groups are labelled
#'     by their values (`"41+42+43"`).
#'   * `"all"` --- every distinct value becomes its own class.
#'
#'   A value may belong to only one class. With more than one class, patch IDs
#'   are unique across classes and `metrics$class` records which class each patch
#'   belongs to.
#' @param directions Connectivity: `8` (default, queen --- diagonal neighbours
#'   are connected) or `4` (rook).
#' @param mask Optional study-area raster, aligned with `x`. Cells that are zero
#'   or `NA` in the mask are outside the analysis domain. Must match `x`'s
#'   geometry exactly; diamondback will not resample silently.
#' @param na How to treat `NA` cells inside the mask. `"outside"` (default)
#'   excludes them from the analysis domain and treats their boundary with a
#'   patch as domain edge. `"background"` treats them as valid background, which
#'   makes them habitat edge.
#' @param crop Crop to the mask's bounding box before labelling, which can cut
#'   memory use substantially when the study area is a small part of the raster.
#'   Defaults to `TRUE` when `mask` is supplied.
#' @param output Optional path to write the labelled raster to. When given, the
#'   raster is written and the result stores the path rather than holding cells
#'   in memory.
#' @param overwrite Overwrite `output` if it exists.
#' @param fingerprint How thoroughly the source is identified for caching and
#'   provenance. `"auto"` (default) hashes file contents under 200 MB and falls
#'   back to size and modification time above that. `"full"` always hashes the
#'   contents, however large: slower, but immune to a file being rewritten at the
#'   same size with its timestamp restored. `"fast"` never hashes. The mode is
#'   recorded and forms part of the cache key, so a result cached under a weaker
#'   fingerprint is never reused for a request that asked for a stronger one.
#' @param max_memory_frac Largest fraction of this machine's **physical** RAM the
#'   arrays may occupy. Above it the operation stops before allocating, because
#'   the job cannot fit however the machine is used. Below it, diamondback will
#'   run and merely warn if memory is currently tight -- an operating system
#'   pages and reclaims, so a job larger than the free memory of the moment is
#'   usually slow rather than impossible, and refusing it would be a false alarm.
#' @param memory_limit Explicit ceiling in bytes, used instead of the detected
#'   physical RAM. For the rare case where you know the budget better than the
#'   machine does, such as a container with a cgroup limit.
#' @param validate Run [validate_patch_result()] afterwards. Costs a second pass;
#'   worth it while setting a workflow up.
#' @param quiet Suppress progress reporting.
#'
#' @return A [patch_result] object. `$patches` gives the labelled `SpatRaster`,
#'   `$metrics` a data frame with one row per patch (ID, class and cell count
#'   only --- use [patch_metrics()] for geometry), and `$metadata` everything
#'   needed to reproduce the run.
#'
#' @seealso [patch_metrics()] for geometry, [analyze_patches()] to do both,
#'   [compare_patches()] for change over time.
#' @export
#' @examplesIf diamondback_ready()
#' # The two blocks touch only at a corner.
#' m <- matrix(c(1, 1, 0, 0,
#'               1, 1, 0, 0,
#'               0, 0, 1, 1,
#'               0, 0, 1, 1), nrow = 4, byrow = TRUE)
#'
#' # Diagonal contact joins them under 8-connectivity, but not under 4.
#' label_patches(m, directions = 8, quiet = TRUE)$metadata$n_patches  # 1
#' label_patches(m, directions = 4, quiet = TRUE)$metadata$n_patches  # 2
#'
#' # A vector is ONE class made of several values; a list is several classes.
#' v <- matrix(c(41, 42, 21, 43), nrow = 2)
#' label_patches(v, class = c(41, 42, 43), quiet = TRUE)$metadata$n_patches      # 1
#' label_patches(v, class = list(41, 42, 43), quiet = TRUE)$metadata$n_patches   # 3
label_patches <- function(x,
                          class = NULL,
                          directions = 8,
                          mask = NULL,
                          na = c("outside", "background"),
                          crop = NULL,
                          output = NULL,
                          overwrite = FALSE,
                          fingerprint = c("auto", "full", "fast"),
                          max_memory_frac = 0.9,
                          memory_limit = NULL,
                          validate = FALSE,
                          quiet = FALSE) {
  na <- match.arg(na)
  fingerprint <- match.arg(fingerprint)
  t_start <- Sys.time()

  if (!is.numeric(directions) || length(directions) != 1L || !directions %in% c(4, 8)) {
    cli::cli_abort(c(
      "{.arg directions} must be {.val 4} or {.val 8}.",
      "x" = "Got {.val {directions}}."
    ), call = NULL)
  }

  src <- db_source_info(x, fingerprint)
  r <- db_as_rast(x)
  mask_r <- db_check_mask(mask, r)
  mask_src <- if (is.null(mask)) NULL else db_source_info(mask, fingerprint)

  if (is.null(crop)) crop <- !is.null(mask_r)
  if (isTRUE(crop) && !is.null(mask_r)) {
    cr <- db_crop_to_mask(r, mask_r, quiet = quiet)
    r <- cr$x
    mask_r <- cr$mask
  }

  class <- db_check_class(class, r)
  if (identical(class, "all")) {
    class <- db_discover_classes(r, quiet = quiet)
    if (!quiet) {
      cli::cli_alert_info("Found {length(class$labels)} class{?es}: {.val {class$labels}}.")
    }
  }
  n_classes <- if (is.null(class)) 1L else length(class$groups)

  geom <- db_geometry(r)
  # More than 253 classes needs a uint16 code array, which is one more byte per
  # cell; the estimate has to know that before anything is allocated.
  code_bytes <- if (n_classes > 253L) 2 else 1
  db_check_memory(geom$ncell, "label", max_memory_frac, memory_limit,
                  quiet = quiet, code_bytes = code_bytes)

  if (!is.null(output)) db_check_output(output, overwrite)

  # ---- ingest ----
  code <- db_build_code(r, mask_r, class, quiet = quiet)
  counts <- db_code_counts(code)

  if (!quiet) {
    cli::cli_alert_info(
      "Cells: {.val {format(counts$foreground, big.mark = ',')}} foreground, \\
       {.val {format(counts$background, big.mark = ',')}} background, \\
       {.val {format(counts$missing, big.mark = ',')}} missing, \\
       {.val {format(counts$outside, big.mark = ',')}} outside domain."
    )
  }

  if (counts$foreground == 0) {
    cli::cli_warn(c(
      "No foreground cells found; the result has zero patches.",
      "i" = if (is.null(class)) {
        "{.arg x} was treated as binary, so every valid cell was zero."
      } else {
        "No cell matched {.arg class}: {.val {unlist(class$groups)}}."
      }
    ))
  }

  # ---- label ----
  py <- db_py()
  use_int64 <- geom$ncell > .Machine$integer.max

  if (!quiet) {
    cli::cli_alert_info("Labelling ({directions}-neighbour connectivity) ...")
  }
  t_label <- Sys.time()

  lab_info <- db_label_classes(py, code, n_classes, directions, use_int64, geom, quiet = quiet)
  labels <- lab_info$labels
  n_total <- lab_info$n_total
  patch_class <- lab_info$patch_class

  db_check_label_range(n_total)

  if (!quiet) {
    cli::cli_alert_success(
      "Found {.val {format(n_total, big.mark = ',')}} \\
       {cli::qty(n_total)}patch{?es} in {db_elapsed(t_label)}."
    )
  }

  # Cell counts fall out of labelling for free and are needed for validation,
  # so they live here rather than in patch_metrics().
  stats <- py_try(py$patch_stats(labels, as.integer(n_total)), "computing patch cell counts")
  cells <- py_num(stats, "count")[-1]

  metrics <- data.frame(
    patch_id = seq_len(n_total),
    cells = cells
  )
  if (n_classes > 1L) {
    # patch_class holds class *indices*; report the label of each class.
    metrics <- data.frame(
      patch_id = metrics$patch_id,
      class = class$labels[patch_class],
      cells = metrics$cells
    )
  }

  # ---- output raster ----
  patches_rast <- db_labels_to_rast(
    py, labels, code, r, output = output, overwrite = overwrite,
    na_background = identical(na, "background"), quiet = quiet
  )

  meta <- db_metadata(
    source = src, mask_source = mask_src, geometry = geom,
    class = class, directions = directions, na = na, crop = crop,
    counts = counts, n_patches = n_total, elapsed = db_elapsed_num(t_start)
  )

  res <- new_patch_result(
    patches = patches_rast$obj,
    patches_path = patches_rast$path,
    metrics = metrics,
    metadata = meta,
    arrays = list(labels = labels, code = code, n = n_total,
                  classes = if (is.null(class)) NULL else class$labels,
                  groups = if (is.null(class)) NULL else class$groups,
                  na_background = identical(na, "background"))
  )

  if (isTRUE(validate)) {
    v <- validate_patch_result(res, quiet = quiet)
    res$metadata$validated <- v$ok
  }

  if (!quiet) cli::cli_alert_success("Done in {db_elapsed(t_start)}.")
  res
}

# ---------------------------------------------------------------------------
# internals
# ---------------------------------------------------------------------------

#' Build the uint8 cell-state code array, one row block at a time
#'
#' The point of the block loop is that the raster is never an R numeric matrix:
#' cells go straight from terra into a preallocated NumPy uint8 array. For a
#' 1e9-cell raster that is 1 GB instead of 8 GB.
#' @noRd
db_build_code <- function(r, mask_r, class, quiet = FALSE) {
  py <- db_py()
  np <- reticulate::import("numpy", convert = FALSE)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  n_classes <- if (is.null(class)) 1L else length(class$groups)

  # uint8 up to 253 classes, uint16 beyond. Chosen once for the whole run from
  # the class count, not per block from whatever values a block contains.
  code <- py_try(np$zeros(reticulate::tuple(as.integer(nr), as.integer(nc)),
                          dtype = py$code_dtype(as.integer(n_classes))),
                 "allocating the cell-state array")

  blocks <- db_row_blocks(nr, nc)
  terra::readStart(r)
  on.exit(terra::readStop(r), add = TRUE)
  if (!is.null(mask_r)) {
    terra::readStart(mask_r)
    on.exit(terra::readStop(mask_r), add = TRUE)
  }

  pb <- db_progress("Reading raster", length(blocks), quiet = quiet)
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    v <- terra::readValues(r, row = b$row, nrows = b$nrows, col = 1, ncols = nc, mat = FALSE)
    mv <- if (is.null(mask_r)) NULL else {
      terra::readValues(mask_r, row = b$row, nrows = b$nrows, col = 1, ncols = nc, mat = FALSE)
    }
    blk <- py_try(
      py$code_block(as.numeric(v), as.integer(b$nrows), as.integer(nc),
                    mask = if (is.null(mv)) NULL else as.numeric(mv),
                    # unname(): reticulate turns a *named* R list into a Python
                    # dict, and enumerate() over a dict yields its keys. The
                    # labels stay on the R side; Python only needs the order.
                    class_groups = if (is.null(class)) NULL else unname(class$groups),
                    n_classes = as.integer(n_classes)),
      "converting raster values to cell states"
    )
    py_try(py$set_rows(code, as.integer(b$row - 1L),
                       as.integer(b$row - 1L + b$nrows), blk),
           "filling the cell-state array")
    db_progress_step(pb, i)
  }
  db_progress_done(pb)
  code
}

#' Row blocks sized so one block is a few million cells
#' @noRd
db_row_blocks <- function(nrow, ncol, target_cells = 4e6) {
  step <- max(1L, as.integer(target_cells %/% max(1L, ncol)))
  starts <- seq.int(1L, nrow, by = step)
  lapply(starts, function(s) list(row = s, nrows = min(step, nrow - s + 1L)))
}

#' Cell-state totals from the code array
#' @noRd
db_code_counts <- function(code) {
  py <- db_py()
  cc <- as.numeric(py_r(py_try(py$code_counts(code), "counting cell states")))
  list(
    outside = cc[1],                 # code 0
    missing = cc[2],                 # code 1
    background = cc[3],              # code 2
    foreground = sum(cc[4:length(cc)]),  # codes 3+
    by_class = cc[4:length(cc)]
  )
}

#' Label each class and merge into globally unique IDs
#' @noRd
db_label_classes <- function(py, code, n_classes, directions, use_int64, geom, quiet = FALSE) {
  if (n_classes == 1L) {
    out <- py_try(
      py$label_array(code, 0L, as.integer(directions), use_int64),
      "connected-component labelling"
    )
    # The label array stays in Python; only the count crosses back.
    return(list(labels = py_get(out, "labels"),
                n_total = as.integer(py_r(py_get(out, "n"))),
                patch_class = NULL))
  }

  labs <- vector("list", n_classes)
  ns <- integer(n_classes)
  offsets <- integer(n_classes)
  total <- 0L
  pb <- db_progress("Labelling classes", n_classes, quiet = quiet)
  for (k in seq_len(n_classes)) {
    out <- py_try(
      py$label_array(code, as.integer(k - 1L), as.integer(directions), use_int64),
      sprintf("connected-component labelling of class %d", k)
    )
    labs[[k]] <- py_get(out, "labels")
    ns[k] <- as.integer(py_r(py_get(out, "n")))
    offsets[k] <- total
    total <- total + ns[k]
    db_progress_step(pb, k)
  }
  db_progress_done(pb)

  combined <- py_try(
    py$combine_labels(labs, as.integer(offsets),
                      reticulate::tuple(as.integer(geom$nrow), as.integer(geom$ncol)),
                      use_int64),
    "merging per-class labels"
  )
  list(labels = combined, n_total = total,
       patch_class = rep(seq_len(n_classes), times = ns))
}

#' Write or materialise the labelled raster
#'
#' Returns either an in-memory SpatRaster or a path, never both, so nothing
#' large is duplicated.
#' @noRd
db_labels_to_rast <- function(py, labels, code, template, output, overwrite,
                              na_background, quiet = FALSE) {
  nr <- terra::nrow(template); nc <- terra::ncol(template)
  out_r <- terra::rast(template)
  terra::varnames(out_r) <- "patch_id"
  names(out_r) <- "patch_id"

  blocks <- db_row_blocks(nr, nc)

  if (is.null(output)) {
    vals <- integer(nr * nc)
    for (b in blocks) {
      rows <- py_try(
        py$output_rows(labels, code, as.integer(b$row - 1L),
                       as.integer(b$row - 1L + b$nrows), na_background),
        "extracting labelled cells"
      )
      idx <- ((b$row - 1L) * nc + 1L):((b$row - 1L + b$nrows) * nc)
      vals[idx] <- db_int_rows(rows)
    }
    terra::values(out_r) <- vals
    return(list(obj = out_r, path = NULL))
  }

  if (!quiet) cli::cli_alert_info("Writing labelled raster to {.path {output}} ...")
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)

  # Write to a sibling temporary file and rename on success, so that a failure
  # part-way through leaves no half-written raster at the user's path -- and, if
  # they were overwriting, leaves the previous file intact rather than truncated.
  # The rename is atomic on the same filesystem, which is why it is a sibling
  # rather than something under tempdir().
  # The marker goes *before* the extension: GDAL infers its driver from the
  # extension, so "labels.tif.part" is a file it cannot write.
  ext <- tools::file_ext(output)
  tmp <- if (nzchar(ext)) {
    paste0(tools::file_path_sans_ext(output), ".diamondback-part.", ext)
  } else {
    paste0(output, ".diamondback-part")
  }
  db_register_file(tmp, kind = "temp")
  unlink(tmp)

  ok <- FALSE
  on.exit({
    if (!ok) {
      unlink(tmp)
      if (!quiet) {
        cli::cli_alert_info("Cleaned up the partial write; {.path {output}} was not modified.")
      }
    }
  }, add = TRUE)

  terra::writeStart(out_r, filename = tmp, overwrite = TRUE,
                    datatype = "INT4S", NAflag = DB_INT32_NA)
  pb <- db_progress("Writing raster", length(blocks), quiet = quiet)
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    rows <- py_try(
      py$output_rows(labels, code, as.integer(b$row - 1L),
                     as.integer(b$row - 1L + b$nrows), na_background),
      "extracting labelled cells"
    )
    terra::writeValues(out_r, db_int_rows(rows), b$row, b$nrows)
    db_progress_step(pb, i)
  }
  terra::writeStop(out_r)
  db_progress_done(pb)

  if (!file.rename(tmp, output)) {
    # Cross-device or a locked target: fall back to a copy, then drop the part.
    if (!file.copy(tmp, output, overwrite = TRUE)) {
      cli::cli_abort(c(
        "Could not move the finished raster into place.",
        "x" = "From {.path {tmp}} to {.path {output}}.",
        "i" = "The completed raster is still at the temporary path."
      ), call = NULL)
    }
    unlink(tmp)
  }
  ok <- TRUE
  .db_files$created <- .db_files$created[.db_files$created != tmp]
  db_register_file(output, kind = "output")

  list(obj = NULL, path = output)
}

#' Crop a raster and mask to the mask's non-zero extent
#' @noRd
db_crop_to_mask <- function(r, mask_r, quiet = FALSE) {
  before <- terra::ncell(r)
  e <- tryCatch({
    mm <- terra::trim(terra::subst(mask_r, 0, NA))
    terra::ext(mm)
  }, error = function(err) NULL)

  if (is.null(e)) return(list(x = r, mask = mask_r))

  e <- terra::intersect(e, terra::ext(r))
  if (is.null(e)) return(list(x = r, mask = mask_r))

  r2 <- terra::crop(r, e, snap = "out")
  m2 <- terra::crop(mask_r, e, snap = "out")
  after <- terra::ncell(r2)

  if (!quiet && after < before) {
    cli::cli_alert_info(
      "Cropped to the mask: {.val {format(before, big.mark = ',')}} -> \\
       {.val {format(after, big.mark = ',')}} cells \\
       ({round(100 * (1 - after / before))}% fewer)."
    )
  }
  list(x = r2, mask = m2)
}

# The largest patch ID that can survive the whole round trip. INT_MIN is R's
# NA_integer_ and is used as the raster's no-data value, so the usable range
# stops at INT_MAX. This is the binding constraint, not SciPy's: scipy can label
# into int64, but an R integer vector, a GeoTIFF INT4S band and terra's NA all
# top out here. Rather than let a 2-billion-patch raster overflow into garbage
# or silent NAs, stop and say so.
MAX_PATCH_ID <- .Machine$integer.max   # 2147483647

#' Refuse to produce patch IDs that cannot be represented
#' @noRd
db_check_label_range <- function(n_total) {
  if (n_total <= MAX_PATCH_ID) return(invisible(TRUE))
  cli::cli_abort(c(
    "Found {format(n_total, big.mark = ',')} patches, which is more than can be \\
     represented in a labelled raster.",
    "x" = "The limit is {format(MAX_PATCH_ID, big.mark = ',')}: patch IDs travel \\
           through R integers and a 32-bit GeoTIFF band, and the value below that \\
           is reserved for {.val NA}.",
    "i" = "A raster this fragmented is usually a sign of noise; try smoothing, \\
           a minimum-mapping-unit filter, or 8-neighbour connectivity.",
    "i" = "Otherwise, split the raster into tiles and label them separately."
  ), class = "diamondback_overflow_error", call = NULL)
}

db_check_output <- function(output, overwrite) {
  if (!is.character(output) || length(output) != 1L) {
    cli::cli_abort("{.arg output} must be a single file path.", call = NULL)
  }
  if (file.exists(output) && !isTRUE(overwrite)) {
    cli::cli_abort(c(
      "{.path {output}} already exists.",
      "i" = "Use {.code overwrite = TRUE} to replace it."
    ), call = NULL)
  }
  invisible(TRUE)
}
