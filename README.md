# diamondback

<!-- badges: start -->
[![R-CMD-check](https://github.com/GatesDupont/diamondback/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GatesDupont/diamondback/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Fast, reliable patch analysis for very large categorical and binary rasters.

`terra::patches()` becomes impractical on rasters with hundreds of millions of
cells. diamondback keeps R in charge of everything that needs to know what a CRS
is — reading, cropping, masking, alignment, metadata — and hands the array
kernels to NumPy and SciPy through reticulate, where connected-component
labelling is orders of magnitude faster.

It exists because the R–Python pipeline for this keeps getting rebuilt one script
at a time, and the hard part was never `scipy.ndimage.label()`. It was the
plumbing around it: `NA` handling, masks, alignment checks, caching, temp files,
and the bookkeeping for patch lineage through time.

## Installation

```r
# install.packages("remotes")
remotes::install_github("GatesDupont/diamondback")
```

Python, NumPy and SciPy are resolved automatically on first use via
`reticulate::py_require()` — nothing is installed into your existing
environments, and nothing happens at `library()` time. To check:

```r
library(diamondback)
diamondback_check()
```

```
-- diamondback environment ----------------------------------------
R: 4.4.1
terra: 1.8.70
reticulate: 1.45.0
Python: 3.12
NumPy: 2.5.1
SciPy: 1.18.0
Available RAM: 21.4 GB

-- End-to-end test --
v Labelling a small test raster succeeded. diamondback is ready.
```

To use your own environment instead, set `RETICULATE_PYTHON` before loading the
package and diamondback will step aside.

## Quick start

```r
library(diamondback)
library(terra)

forest <- rast("forest_1985.tif")
study  <- rast("study_area.tif")

p1985 <- analyze_patches(
  forest,
  class      = 1,
  mask       = study,
  directions = 8,
  output_dir = "derived/1985"   # writes results here, and caches them
)

p1985
#> -- <patch_result> -------------------------------------------------
#> 84,213 patches from 412,000,000 cells (20000 x 20600)
#> source: forest_1985.tif
#> class: 1
#> connectivity: 8-neighbour
#> NA policy: outside
#> ...
```

Three things to reach for:

```r
p1985$patches     # the labelled SpatRaster (from disk or memory, transparently)
p1985$metrics     # one row per patch: area, perimeter, edges, bbox, centroid
p1985$metadata    # everything needed to reproduce the run
```

Then compare two times:

```r
p2024  <- analyze_patches(rast("forest_2024.tif"), class = 1, mask = study,
                          output_dir = "derived/2024")

change <- compare_patches(p1985, p2024)

change$overlaps   # every overlapping pair, with proportions and ranks
change$events     # one row per patch: split / merger / persistence / ...
change$lineages   # one row per connected group of related patches
```

Re-running `analyze_patches()` with the same inputs reads the cached result
instead of relabelling. Change the class, the mask, the connectivity or the
source file, and it recomputes and tells you which one changed.

## Is it actually faster?

Measured on an M-series Mac at 40% forest cover with 8-connectivity
(`inst/benchmarks/benchmark-terra.R`), checking at every size that the two
agree on the patch count:

| raster | cells | diamondback | `terra::patches()` | speedup |
| --- | --- | --- | --- | --- |
| 500×500 | 250,000 | 0.18 s | 0.55 s | 3× |
| 1000×1000 | 1,000,000 | 0.51 s | 8.34 s | ~16× |
| 2000×2000 | 4,000,000 | 0.96 s | **did not finish in 20 min** | >1000× |

`terra::patches()` degrades superlinearly as the patch count grows. diamondback
stays roughly linear in cells:

| raster | cells | label | metrics | core area | patches found |
| --- | --- | --- | --- | --- | --- |
| 4000×4000 | 16,000,000 | 2.7 s | 1.9 s | 2.3 s | 360,866 |
| 8000×8000 | 64,000,000 | 12.4 s | 8.8 s | 11.2 s | 1,444,076 |

Don't read the speedup as a fixed number — it grows with raster size and patch
count, which is exactly why this package exists. On a small raster, `terra` is
perfectly fine and has no Python dependency.

## The part that matters most: `NA`

The most common way these analyses go quietly wrong is an `NA` becoming
foreground or background by accident. diamondback keeps four states apart and
never merges them silently:

| state | meaning | in the labelled raster |
| --- | --- | --- |
| **patch** | a foreground cell | `1..N` |
| **background** | valid, in-domain, not the target class | `0` |
| **outside domain** | excluded by your `mask` | `NA` |
| **missing** | `NA` in the source, inside the mask | `NA` |

`NA` is never foreground. It is not background either unless you say so with
`na = "background"`. Exact counts for all four are in
`result$metadata$cells`, and `patch_domain()` returns them as a categorical
raster.

## Edges: habitat versus study area

`perimeter` splits into two parts that always sum back to it:

- `edge_valid` — boundary against a real, in-domain cell. The patch genuinely
  stops here.
- `edge_domain` — boundary against a masked cell, a missing cell, or the grid
  border. The patch may well continue; you just cannot see it.

`touches_domain_edge` flags patches whose area and perimeter are therefore
**lower bounds** — usually the ones to drop from a size distribution.

Perimeter counts shared cell edges only. Corner contact contributes nothing,
because two cells touching at a corner share no physical boundary. Non-square
cells are handled exactly: north–south edges use `yres`, east–west use `xres`.

## Core area

```r
p <- patch_core_area(p1985, depth = 300, units = "m")
p$metrics[, c("patch_id", "cells", "core_cells", "core_fraction")]
```

One exact Euclidean distance transform, any depth, no repeated erosion.
Distance is measured centre-to-centre, so a cell adjacent to non-habitat is one
cell width from the edge, not zero — see `?patch_core_area` for the precise rule
and for `edge = "background"`, which spares patches clipped by the study
boundary.

## Time series

```r
series <- track_patch_series(
  c("forest_1985.tif", "forest_2000.tif", "forest_2024.tif"),
  time       = c(1985, 2000, 2024),
  class      = 1,
  mask       = study,
  output_dir = "derived"
)

series$summary       # patch counts and event tallies per transition
series$transitions   # every patch, every transition, stacked
```

Each step is cached, so an interrupted run resumes instead of relabelling.

## Large rasters

Labelling holds the array in RAM — that is `scipy.ndimage.label()`'s contract,
and diamondback is honest about it rather than pretending otherwise. What it does
do:

```r
db_memory_report(rast("big.tif"))
#> Raster: 20000 x 20600 = 412,000,000 cells
#> label_patches()   2.3 GB
#> patch_metrics()   1.9 GB
#> patch_core_area() 5.0 GB
#> Available RAM     21.4 GB
```

Every entry point estimates before it allocates and **errors before**, not
during, an allocation that will not fit. Cells are read blockwise straight into
a `uint8` NumPy array, so the raster is never an R numeric matrix (an 8×
saving). With a mask, the grid is cropped to it first.

diamondback also never calls `terra::tmpFiles(remove = TRUE)`. It tracks the
files it creates and removes only those — deleting files that back someone
else's live `SpatRaster` is exactly the failure this package is meant to end.

## Trust

```r
validate_patch_result(p1985)
#> v All 9 validation checks passed.
```

Independent re-derivations, not assertions restating the code: patch counts
recounted from the raster, foreground recounted from the cell states, sums
compared, mask leakage checked, edge components reconciled against perimeter.

The test suite checks diamondback against `terra::patches()` on small rasters
under both connectivity rules, at several foreground densities, and with `NA`
present — the two must agree exactly on the partition, not just the count.

## API

| function | purpose |
| --- | --- |
| `diamondback_check()` | environment diagnostic, end to end |
| `analyze_patches()` | label + measure, with caching — **start here** |
| `label_patches()` | connected-component labelling |
| `patch_metrics()` | area, perimeter, edges, bbox, centroid |
| `patch_core_area()` | habitat beyond a distance from an edge |
| `compare_patches()` | overlap, lineage and events between two times |
| `track_patch_series()` | the same across a time series |
| `patch_domain()` | the four cell states as a raster |
| `validate_patch_result()` | independent correctness checks |
| `write_patch_result()` / `read_patch_result()` | durable results on disk |
| `db_memory_report()` | what a raster will cost before you commit |
| `db_clean_temp()` | remove only diamondback's own temp files |

## Design

[`DESIGN.md`](DESIGN.md) documents the architecture and the reasoning behind
each decision: why labelling and metrics are separate, how perimeter is defined,
how the mother patch is chosen, what invalidates a cache, and what was
considered and rejected for v1.

`vignette("diamondback")` walks through a complete workflow.

## Limitations

- Labelling requires the array to fit in RAM. Tiled/Dask backends were
  considered and rejected for v1; the kernel is isolated so one can be added
  without changing the API.
- `patch_core_area()` needs a projected raster. Everything else supports lon/lat
  exactly, using per-row geodesic geometry.
- Up to 252 classes in one run.
- Lineage is built from consecutive pairwise comparisons, not a global identity
  model.

## License

MIT © Gates Dupont
