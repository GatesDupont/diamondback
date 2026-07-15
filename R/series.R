# Multi-year tracking. Version 1 chains consecutive comparisons rather than
# attempting a global identity model; the data structures are shaped so that a
# better lineage model can be layered on without changing the API.

#' Track patches across a time series
#'
#' Labels each time step, compares consecutive steps, and stacks the results
#' into one long lineage table. This is deliberately a chain of pairwise
#' comparisons: it does not try to solve persistent patch identity across the
#' whole series, which is a genuinely unsettled problem. What it gives you is
#' every consecutive transition in one tidy table, which is enough to build most
#' lineage analyses and is honest about what it knows.
#'
#' @param x A multi-layer `SpatRaster`, a vector of raster filenames, or a list
#'   of `SpatRaster`s or [patch_result]s.
#' @param time Labels for the time steps: a vector as long as `x`. Defaults to
#'   layer names, file basenames, or `1..n`.
#' @param class,directions,mask,na,crop Passed to [label_patches()]. Note that
#'   `crop` applies identically to every step, which is what keeps the series
#'   comparable: cropping each year to its own extent would misalign them.
#' @param min_overlap_cells,min_overlap_prop Passed to [compare_patches()].
#' @param output_dir Optional directory. When given, each time step is cached
#'   there (see [analyze_patches()]), so an interrupted run resumes rather than
#'   relabelling everything.
#' @param metrics Also compute [patch_metrics()] for each time step.
#' @param quiet Suppress progress reporting.
#'
#' @return A `patch_series` object with:
#'   * `$results` --- named list of [patch_result]s, one per time step.
#'   * `$transitions` --- one row per patch per transition, stacked `$events`
#'     tables with `time_from` and `time_to` columns added.
#'   * `$overlaps` --- stacked `$overlaps` tables with the same two columns.
#'   * `$summary` --- one row per transition: patch counts, cell totals, and
#'     event tallies.
#'
#' @seealso [compare_patches()]
#' @export
#' @examplesIf diamondback_ready()
#' # Three steps: a patch grows, then splits.
#' m1 <- matrix(0, 6, 6); m1[2:3, 2:3] <- 1
#' m2 <- matrix(0, 6, 6); m2[2:5, 2:5] <- 1
#' m3 <- matrix(0, 6, 6); m3[2:3, 2:5] <- 1; m3[5, 2:5] <- 1
#' s <- track_patch_series(list(m1, m2, m3), time = c(2000, 2010, 2020), quiet = TRUE)
#' s$summary
track_patch_series <- function(x,
                               time = NULL,
                               class = NULL,
                               directions = 8,
                               mask = NULL,
                               na = c("outside", "background"),
                               crop = NULL,
                               min_overlap_cells = 1L,
                               min_overlap_prop = 0,
                               output_dir = NULL,
                               metrics = TRUE,
                               quiet = FALSE) {
  na <- match.arg(na)
  t0 <- Sys.time()

  steps <- db_series_steps(x)
  time <- db_series_time(time, steps, x)
  n <- length(steps)

  if (n < 2L) {
    cli::cli_abort(c(
      "A series needs at least two time steps.",
      "x" = "Got {n}.",
      "i" = "For a single raster use {.fn analyze_patches}."
    ), call = NULL)
  }

  results <- vector("list", n)
  names(results) <- as.character(time)

  for (i in seq_len(n)) {
    if (!quiet) cli::cli_h2("Time step {i}/{n}: {time[i]}")
    if (inherits(steps[[i]], "patch_result")) {
      results[[i]] <- steps[[i]]
      next
    }
    od <- if (is.null(output_dir)) NULL else file.path(output_dir, as.character(time[i]))
    results[[i]] <- analyze_patches(
      steps[[i]], class = class, directions = directions, mask = mask, na = na,
      crop = crop, output_dir = od, metrics = metrics, quiet = quiet
    )
  }

  transitions <- list()
  overlaps <- list()
  summ <- list()

  for (i in seq_len(n - 1L)) {
    if (!quiet) cli::cli_h2("Transition {i}/{n - 1L}: {time[i]} -> {time[i + 1L]}")
    cmp <- compare_patches(
      results[[i]], results[[i + 1L]],
      min_overlap_cells = min_overlap_cells,
      min_overlap_prop = min_overlap_prop,
      quiet = quiet
    )

    ev <- cmp$events
    ev$time_from <- time[i]
    ev$time_to <- time[i + 1L]
    transitions[[i]] <- ev

    ov <- cmp$overlaps
    if (nrow(ov)) {
      ov$time_from <- time[i]
      ov$time_to <- time[i + 1L]
    }
    overlaps[[i]] <- ov

    tb <- table(factor(cmp$lineages$event,
                       levels = c("persistence", "split", "merger", "complex",
                                  "appearance", "disappearance")))
    summ[[i]] <- data.frame(
      time_from = time[i],
      time_to = time[i + 1L],
      n_patches_from = cmp$metadata$n_patches_t1,
      n_patches_to = cmp$metadata$n_patches_t2,
      cells_from = sum(ev$cells[ev$time == "t1"]),
      cells_to = sum(ev$cells[ev$time == "t2"]),
      persistence = as.integer(tb[["persistence"]]),
      split = as.integer(tb[["split"]]),
      merger = as.integer(tb[["merger"]]),
      complex = as.integer(tb[["complex"]]),
      appearance = as.integer(tb[["appearance"]]),
      disappearance = as.integer(tb[["disappearance"]])
    )
  }

  res <- structure(
    list(
      results = results,
      transitions = do.call(rbind, transitions),
      overlaps = do.call(rbind, overlaps),
      summary = do.call(rbind, summ),
      metadata = list(
        created = as.character(Sys.time()),
        package_version = as.character(utils::packageVersion("diamondback")),
        elapsed_secs = db_elapsed_num(t0),
        time = time,
        n_steps = n,
        directions = directions,
        class = if (is.null(class)) NA_real_ else class,
        min_overlap_cells = min_overlap_cells,
        min_overlap_prop = min_overlap_prop
      )
    ),
    class = "patch_series"
  )
  if (!quiet) cli::cli_alert_success("Tracked {n} time step{?s} in {db_elapsed(t0)}.")
  res
}

db_series_steps <- function(x) {
  if (inherits(x, "SpatRaster")) {
    if (terra::nlyr(x) < 2L) {
      cli::cli_abort("{.arg x} is a single-layer raster; a series needs at least two.", call = NULL)
    }
    return(lapply(seq_len(terra::nlyr(x)), function(i) x[[i]]))
  }
  if (is.character(x)) {
    missing <- x[!file.exists(x)]
    if (length(missing)) {
      cli::cli_abort(c(
        "{length(missing)} file{?s} in {.arg x} do{?es/} not exist.",
        "x" = "{.path {utils::head(missing, 5)}}"
      ), call = NULL)
    }
    return(as.list(x))
  }
  if (is.list(x)) return(x)
  cli::cli_abort(c(
    "{.arg x} must be a multi-layer {.cls SpatRaster}, a vector of filenames, or a list.",
    "x" = "Got {.cls {class(x)[1]}}."
  ), call = NULL)
}

db_series_time <- function(time, steps, x) {
  n <- length(steps)
  if (!is.null(time)) {
    if (length(time) != n) {
      cli::cli_abort(c(
        "{.arg time} must have one entry per time step.",
        "x" = "Got {length(time)} for {n} step{?s}."
      ), call = NULL)
    }
    if (anyDuplicated(time)) {
      cli::cli_abort("{.arg time} has duplicate values; time steps must be distinct.", call = NULL)
    }
    return(time)
  }
  if (inherits(x, "SpatRaster")) return(names(x))
  if (is.character(x)) return(tools::file_path_sans_ext(basename(x)))
  if (!is.null(names(steps))) return(names(steps))
  seq_len(n)
}

#' @export
print.patch_series <- function(x, ...) {
  m <- .subset2(x, "metadata")
  cli::cli_h1("<patch_series>")
  cli::cli_dl(c(
    "time steps" = "{m$n_steps} ({paste(utils::head(m$time, 6), collapse = ', ')}\\
                     {if (m$n_steps > 6) ', ...' else ''})",
    "connectivity" = "{m$directions}-neighbour",
    "transitions" = "{nrow(.subset2(x, 'summary'))}"
  ))
  s <- .subset2(x, "summary")
  if (nrow(s)) {
    cli::cli_h2("Transitions")
    print(s[, c("time_from", "time_to", "n_patches_from", "n_patches_to",
                "split", "merger", "appearance", "disappearance")],
          row.names = FALSE)
  }
  cli::cli_text("")
  cli::cli_text("{.strong $results} {.strong $transitions} {.strong $overlaps} {.strong $summary}")
  invisible(x)
}
