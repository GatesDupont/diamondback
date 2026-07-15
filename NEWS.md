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

* **Running an analysis never downloads or installs anything.** diamondback uses
  a Python you already have — `RETICULATE_PYTHON`, an active virtualenv or conda
  environment, or any interpreter on `PATH` that has NumPy and SciPy. If it finds
  none, it stops and explains the options rather than fetching an interpreter.
* `diamondback_install_python()` is the only thing that downloads. It says what
  it will fetch, asks first when interactive, and records consent so it asks
  once. `diamondback_remove_python()` revokes it.
* `diamondback_check()` reports the full stack, names which Python source was
  used, and runs a real labelling operation end to end.
* Most users need configure nothing: any conda environment or virtualenv with
  NumPy and SciPy is found automatically, deterministically, and reported.
* `options(diamondback.python = )` accepts an environment **name**
  (`"geo"`, `"python3"`) as well as a path, resolved against active
  environments, the `PATH`, conda installations and virtualenvs. Names are
  verified to have NumPy and SciPy before use.
* When no suitable Python exists, the error lists the environments that were
  found rather than leaving you to hunt for a path.

## Limits and conventions made explicit

* Up to 65,533 classes per run. The cell-state array widens from `uint8` to
  `uint16` automatically past 253 classes, at one extra byte per cell, which the
  memory estimate accounts for.
* `fingerprint = c("auto", "full", "fast")` controls how strongly a source is
  identified for caching. `"auto"` hashes files under 200 MB and uses size and
  mtime above that; `"full"` always hashes, for pipelines where inputs are
  regenerated and a size-and-mtime match could be a false hit. The mode is part
  of the cache key.
* Rotated rasters and lon/lat extents beyond the poles are rejected rather than
  approximated. The supported geometry is a north-up, axis-aligned, regular grid.
* Patch IDs are deterministic labels, not persistent identities: never join two
  runs on `patch_id`.
* Single-process and single-threaded; the kernels never touch BLAS, so results
  are bit-for-bit reproducible and thread-count environment variables have no
  effect.
* Core-area distance is centre-to-centre, documented and pinned by tests. A
  nominal depth is measured to the nearest non-habitat cell *centre*, half a cell
  further than the physical patch boundary.
