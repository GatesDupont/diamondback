#' @keywords internal
#' @aliases diamondback-package
#'
#' @section Getting started:
#' Check the environment once, then analyse:
#'
#' ```r
#' diamondback_check()
#'
#' library(terra)
#' forest <- rast("forest_1985.tif")
#' study  <- rast("study_area.tif")
#'
#' p1985 <- analyze_patches(forest, class = 1, mask = study, output_dir = "derived/1985")
#' p2024 <- analyze_patches(rast("forest_2024.tif"), class = 1, mask = study,
#'                          output_dir = "derived/2024")
#'
#' change <- compare_patches(p1985, p2024)
#' ```
#'
#' @section Core functions:
#' * [analyze_patches()] --- label patches and measure them, with caching. Start here.
#' * [label_patches()] --- connected-component labelling on its own.
#' * [patch_metrics()] --- area, perimeter, edge composition, bbox, centroid.
#' * [patch_core_area()] --- habitat beyond a given distance from an edge.
#' * [compare_patches()] --- overlap, lineage and events between two times.
#' * [track_patch_series()] --- the same across a time series.
#' * [validate_patch_result()] --- independent correctness checks.
#' * [write_patch_result()] / [read_patch_result()] --- durable results on disk.
#' * [diamondback_check()] --- environment diagnostic.
#'
#' @section How NA is treated:
#' diamondback keeps four cell states apart --- patch, valid background, outside
#' the analysis domain, and missing --- and `NA` is never silently folded into
#' any of them. See [label_patches()] for the full account. This is the single
#' most important thing to understand before trusting the numbers.
#'
#' @section Python:
#' The array kernels run in NumPy and SciPy via reticulate. Nothing is installed
#' when the package loads; the backend is resolved the first time a computation
#' needs it, using `reticulate::py_require()`. Run [diamondback_check()] if
#' anything looks wrong.
#'
#' @section Limitations in this version:
#' Labelling holds the array in RAM, which is what `scipy.ndimage.label()`
#' requires. [db_memory_report()] estimates the cost of a raster before you
#' commit to it, and every entry point errors before allocating rather than
#' after. [patch_core_area()] needs a projected raster; everything else supports
#' lon/lat exactly.
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
