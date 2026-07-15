# Independent checks on a finished result. These are cheap relative to
# labelling and catch exactly the class of bug the old pipelines shipped
# silently: an NA that became foreground, a mask that leaked, a table that
# disagrees with its raster.

#' Validate a patch result
#'
#' Re-derives facts about a result and checks them against what the result
#' claims. These are independent checks, not assertions restating the code that
#' produced the numbers: patch counts are recounted from the raster, foreground
#' cells are recounted from the cell states, and the sums are compared.
#'
#' Run this while setting a workflow up, and on the first real run of a long
#' pipeline. `label_patches(validate = TRUE)` wires it in automatically.
#'
#' @section What is checked:
#' * every patch ID from `1..N` is present and consecutive;
#' * the foreground cell count equals the sum of patch sizes;
#' * the number of patches in `$metrics` equals the number in the raster;
#' * no foreground cell lies outside the mask;
#' * no patch has zero cells;
#' * patch IDs fit the integer type used;
#' * `edge_valid + edge_missing + edge_outside` equals `perimeter`;
#' * core area never exceeds total area.
#'
#' @param x A [patch_result].
#' @param error Raise an error on failure instead of returning the report.
#' @param quiet Suppress the report.
#' @return A list with `ok` (logical) and `checks` (a data frame with one row per
#'   check: `check`, `passed`, `detail`), invisibly.
#' @export
#' @examplesIf diamondback_ready()
#' m <- matrix(c(1, 1, 0, NA, 1, 0, 0, 1, 1), nrow = 3)
#' v <- validate_patch_result(analyze_patches(m, quiet = TRUE))
#' v$ok
validate_patch_result <- function(x, error = FALSE, quiet = FALSE) {
  stopifnot(inherits(x, "patch_result"))
  meta <- .subset2(x, "metadata")
  met <- .subset2(x, "metrics")
  checks <- list()

  add <- function(name, passed, detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, passed = isTRUE(passed), detail = detail
    )
  }

  st <- tryCatch(db_result_state(x, quiet = TRUE), error = function(e) NULL)
  if (is.null(st)) {
    add("arrays available", FALSE, "could not rebuild arrays from the result")
    return(db_validate_finish(checks, error, quiet))
  }

  py <- db_py()
  n <- st$n

  # --- patch count agreement between raster and table ---
  add("metrics rows match patch count",
      nrow(met) == n,
      sprintf("%d rows, %d patches in raster", nrow(met), n))

  # --- recount cells independently ---
  stats <- py_try(py$patch_stats(st$labels, as.integer(n)), "revalidating patch counts")
  cells <- py_num(stats, "count")[-1]
  counts <- db_code_counts(st$code)

  add("patch IDs are consecutive from 1",
      n == 0L || (all(cells > 0)),
      if (n > 0 && any(cells == 0)) {
        sprintf("%d ID(s) with zero cells, e.g. %s",
                sum(cells == 0), paste(utils::head(which(cells == 0), 5), collapse = ", "))
      } else "")

  add("foreground cells equal summed patch sizes",
      isTRUE(all.equal(sum(cells), counts$foreground)),
      sprintf("%s foreground vs %s summed",
              format(counts$foreground, big.mark = ","), format(sum(cells), big.mark = ",")))

  if ("cells" %in% names(met) && nrow(met) == n && !anyNA(met$cells)) {
    add("table cell counts match raster",
        isTRUE(all.equal(as.numeric(met$cells), cells)),
        if (!isTRUE(all.equal(as.numeric(met$cells), cells))) {
          sprintf("%d patch(es) disagree", sum(as.numeric(met$cells) != cells))
        } else "")
  }

  # --- no foreground outside the mask ---
  # By construction the code array cannot hold both states, so this recomputes
  # it from the array rather than trusting the construction.
  leak <- tryCatch({
    n_cls <- if (is.null(st$classes)) 1L else length(st$classes)
    sum(vapply(seq_len(n_cls), function(k) {
      as.numeric(py_r(py_try(py$foreground_outside_domain(st$code, as.integer(k - 1L)),
                             "checking for foreground outside the mask")))
    }, numeric(1)))
  }, error = function(e) NA_real_)
  add("no foreground outside the analysis domain",
      !is.na(leak) && leak == 0,
      if (!is.na(leak) && leak > 0) sprintf("%s cell(s) leaked", format(leak, big.mark = ",")) else "")

  # --- integer range ---
  add("patch IDs fit in a 32-bit integer",
      n <= .Machine$integer.max,
      sprintf("%s patches", format(n, big.mark = ",")))

  # --- cell-state totals reconcile ---
  tot <- counts$foreground + counts$background + counts$missing + counts$outside
  add("cell states account for every cell",
      isTRUE(all.equal(tot, meta$geometry$ncell)),
      sprintf("%s counted vs %s in raster",
              format(tot, big.mark = ","), format(meta$geometry$ncell, big.mark = ",")))

  # --- metric internal consistency ---
  if (all(c("perimeter", "edge_valid", "edge_missing", "edge_outside") %in% names(met)) &&
      nrow(met)) {
    d <- abs(met$edge_valid + met$edge_missing + met$edge_outside - met$perimeter)
    tol <- 1e-6 * pmax(1, met$perimeter)
    add("edge components sum to perimeter",
        all(d <= tol),
        if (any(d > tol)) sprintf("%d patch(es) disagree", sum(d > tol)) else "")
  }
  if ("core_cells" %in% names(met) && "cells" %in% names(met) && nrow(met)) {
    add("core cells never exceed patch cells",
        all(met$core_cells <= met$cells),
        if (any(met$core_cells > met$cells)) {
          sprintf("%d patch(es) exceed", sum(met$core_cells > met$cells))
        } else "")
  }
  if ("area_m2" %in% names(met) && nrow(met)) {
    add("all patch areas are positive",
        all(met$area_m2 > 0),
        if (any(met$area_m2 <= 0)) sprintf("%d non-positive", sum(met$area_m2 <= 0)) else "")
  }

  db_validate_finish(checks, error, quiet)
}

db_validate_finish <- function(checks, error, quiet) {
  tab <- do.call(rbind, checks)
  ok <- all(tab$passed)

  if (!quiet) {
    if (ok) {
      cli::cli_alert_success("All {nrow(tab)} validation check{?s} passed.")
    } else {
      failed <- tab[!tab$passed, ]
      cli::cli_alert_danger("{nrow(failed)} of {nrow(tab)} validation check{?s} failed:")
      for (i in seq_len(nrow(failed))) {
        cli::cli_bullets(c("x" = "{failed$check[i]}{if (nzchar(failed$detail[i])) paste0(': ', failed$detail[i]) else ''}"))
      }
    }
  }

  if (!ok && isTRUE(error)) {
    failed <- tab[!tab$passed, ]
    cli::cli_abort(c(
      "Patch result failed validation.",
      stats::setNames(
        paste0(failed$check, ifelse(nzchar(failed$detail), paste0(": ", failed$detail), "")),
        rep("x", nrow(failed))
      ),
      "i" = "This should not happen; please report it with a reproducible example."
    ), class = "diamondback_validation_error", call = NULL)
  }

  invisible(list(ok = ok, checks = tab))
}
