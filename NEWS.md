# diamondback 0.1.0

First release.

## Core

* `label_patches()` — connected-component labelling via `scipy.ndimage.label()`,
  for binary rasters, a selected class of a categorical raster, or every class
  separately. 4- and 8-neighbour connectivity, study-area masks, and explicit
  `NA` handling.
* `patch_metrics()` — area, perimeter, edge composition, bounding box and
  centroid. Non-square cells and geographic rasters handled exactly.
* `patch_core_area()` — habitat beyond a given edge depth, via one exact
  Euclidean distance transform.
* `compare_patches()` — overlap, lineage and event classification between two
  times.
* `track_patch_series()` — the same across a time series.
* `analyze_patches()` — label and measure in one call, with caching.

## Correctness

* Four cell states are kept distinct throughout: patch, valid background,
  outside the analysis domain, and missing. `NA` is never silently folded into
  foreground or background. `patch_domain()` returns the distinction as a
  raster.
* `perimeter` splits into `edge_valid` (real habitat edge) and `edge_domain`
  (the study-area boundary or grid border), and `touches_domain_edge` flags
  patches whose metrics are lower bounds.
* `validate_patch_result()` independently re-derives patch counts, foreground
  totals, mask containment and metric consistency.
* Verified against `terra::patches()` on small rasters under both connectivity
  rules, at several foreground densities, and with `NA` present.

## Safety

* Memory is estimated before allocation; operations error *before* an
  allocation that will not fit, and `db_memory_report()` reports the cost of a
  raster without doing any work.
* Raster cells are read blockwise straight into a `uint8` NumPy array, so a
  raster is never materialised as an R numeric matrix.
* `terra::tmpFiles(remove = TRUE)` is never called. Only files diamondback
  created are removed, and `db_clean_temp()` leaves your `output` files alone.
* `compare_patches()` and `mask` refuse to resample silently; misaligned inputs
  are an error naming the mismatched property.

## Environment

* Python, NumPy and SciPy are resolved on first use via
  `reticulate::py_require()`. Nothing is installed into existing environments,
  and nothing happens at `library()` time.
* `diamondback_check()` reports the full stack and runs a real labelling
  operation end to end.
