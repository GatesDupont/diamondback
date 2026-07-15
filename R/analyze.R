#' Label patches and calculate their metrics
#'
#' The function most workflows should call. Combines [label_patches()] and
#' [patch_metrics()], and adds caching: with `output_dir`, the result is written
#' there, and a later call with the same inputs reads it back instead of
#' relabelling.
#'
#' @section Caching:
#' A cached result is reused only when it is provably the same computation. The
#' key covers the source fingerprint (path, size, modification time, and content
#' hash for files under 200 MB), the raster geometry, the class, the mask, the
#' connectivity, the NA policy, and the labelling algorithm version. Change any
#' of them and the work is redone, with a message naming the field that differed.
#'
#' Things that deliberately do **not** invalidate a cache: metric options, units,
#' progress settings, and the package version. None of them can change the
#' labels.
#'
#' A large in-memory raster with no file behind it cannot be fingerprinted, and
#' is therefore never served from cache. Point `x` at a file to get caching.
#'
#' @inheritParams label_patches
#' @param output_dir Directory for the result. When set, the labelled raster,
#'   metrics and metadata are written there, and a matching cached result is
#'   reused. When `NULL` (default) nothing is written and nothing is cached.
#' @param cache Cache behaviour: `TRUE` (default) uses a matching cached result,
#'   `FALSE` recomputes and overwrites, `"read"` uses the cache but does not
#'   write.
#' @param metrics Compute [patch_metrics()]. `FALSE` gives labels and cell
#'   counts only.
#' @param units,edge_depth Passed to [patch_metrics()].
#'
#' @return A [patch_result].
#'
#' @seealso [label_patches()], [patch_metrics()], [compare_patches()]
#' @export
#' @examples
#' m <- matrix(c(1, 1, 0, 0,
#'               1, 1, 0, 1,
#'               0, 0, 0, 1,
#'               0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' res <- analyze_patches(m, quiet = TRUE)
#' res$metrics[, c("patch_id", "cells", "perimeter", "touches_domain_edge")]
analyze_patches <- function(x,
                            class = NULL,
                            directions = 8,
                            mask = NULL,
                            na = c("outside", "background"),
                            crop = NULL,
                            output_dir = NULL,
                            cache = TRUE,
                            metrics = TRUE,
                            units = c("m", "km", "cells"),
                            edge_depth = NULL,
                            overwrite = FALSE,
                            max_memory_frac = 0.6,
                            memory_limit = NULL,
                            validate = FALSE,
                            quiet = FALSE) {
  na <- match.arg(na)
  units <- match.arg(units)
  t0 <- Sys.time()

  if (!is.null(output_dir) && !isFALSE(cache)) {
    req <- db_request_meta(x, class, directions, mask, na, crop)
    hit <- db_cache_lookup(output_dir, req, quiet = quiet)
    if (!is.null(hit)) {
      # A cached run without metrics still needs them if this call wants them.
      if (isTRUE(metrics) && !("perimeter" %in% names(hit$metrics))) {
        hit <- patch_metrics(hit, units = units, edge_depth = edge_depth, quiet = quiet)
      }
      return(hit)
    }
  }

  out_tif <- if (is.null(output_dir)) NULL else file.path(output_dir, RESULT_FILES$labels)

  res <- label_patches(
    x, class = class, directions = directions, mask = mask, na = na, crop = crop,
    output = out_tif, overwrite = overwrite || !is.null(output_dir),
    max_memory_frac = max_memory_frac, memory_limit = memory_limit,
    validate = FALSE, quiet = quiet
  )

  if (isTRUE(metrics)) {
    res <- patch_metrics(res, units = units, edge_depth = edge_depth, quiet = quiet)
  }

  if (isTRUE(validate)) {
    v <- validate_patch_result(res, quiet = quiet)
    res$metadata$validated <- v$ok
  }

  if (!is.null(output_dir) && !identical(cache, "read")) {
    write_patch_result(res, output_dir, overwrite = TRUE, quiet = quiet)
  }

  res$metadata$elapsed_secs <- db_elapsed_num(t0)
  res
}
