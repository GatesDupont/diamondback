# diamondback: design notes

Version 1 design rationale and the decisions behind the API. This is the document
to read before changing the computational core.

## 1. Problem

`terra::patches()` is impractically slow on rasters with hundreds of millions of
cells. `scipy.ndimage.label()` completes the same connected-component labeling in
a fraction of the time on an in-memory array. The recurring cost is not the
algorithm, it is the *plumbing*: NA handling, masks, geometry validation,
temp-file management, caching, and the tabular bookkeeping around lineage.

diamondback keeps R in charge of geometry, IO and validation, and hands only the
array kernels to Python.

## 2. Architecture

```
                    R                                Python
 ┌──────────────────────────────────────┐   ┌──────────────────────────┐
 │ input normalisation  (input.R)       │   │                          │
 │ geometry + alignment validation      │   │ diamondback_core.py      │
 │ memory estimate       (memory.R)     │   │  code_block()            │
 │ blockwise raster read (terra)        │──▶│  label_array()   scipy   │
 │ progress / cli                       │   │  patch_stats()   numpy   │
 │ metric assembly       (metrics.R)    │◀──│  edge_lengths()          │
 │ lineage tables        (compare.R)    │   │  core_counts()   EDT     │
 │ caching + fingerprints(cache.R)      │   │  overlap_counts()        │
 │ S3 results            (result.R)     │   │                          │
 └──────────────────────────────────────┘   └──────────────────────────┘
```

The Python module is shipped in `inst/python/` and loaded with
`reticulate::import_from_path()`. Users never source it.

**Division of labour.** Python does exactly the things that need an array kernel:
connected-component labeling, `bincount`-style reductions, boundary-edge
accumulation, the Euclidean distance transform, and the overlap cross-tabulation.
Everything that needs to know what a CRS is stays in R. The rule of thumb: if the
operation is O(cells) and index-based, it belongs in Python; if it is O(patches)
or needs geometry, it belongs in R.

## 3. The four cell states

The single most common source of silent error in the old pipelines was `NA`
becoming foreground or background by accident. diamondback therefore carries an
explicit **class code array** (`uint8`, one byte per cell) alongside the labels:

| code   | meaning                                       |
|--------|-----------------------------------------------|
| `0`    | outside the analysis domain (mask is FALSE/NA) |
| `1`    | missing / unknown (`NA` in the source, inside the mask) |
| `2`    | valid background (in domain, not the target class) |
| `3+k`  | foreground of class *k*                        |

Precedence when building the codes: outside-domain beats missing, missing beats
class assignment. So a cell that is `NA` in the source *and* outside the mask is
`0`, not `1` — outside the domain we do not care why.

`NA` is **never** foreground. By default it is not background either: it is its
own state, and it is reported separately in `metadata$cells`. Users who genuinely
mean "NA is absence" opt in explicitly with `na = "background"`.

**In the returned raster** the four states collapse to three, because a raster has
only one no-data value:

| labelled raster value | meaning                       |
|-----------------------|-------------------------------|
| `1..N`                | patch ID                       |
| `0`                   | valid background               |
| `NA`                  | outside domain *or* missing    |

The distinction between "outside" and "missing" is preserved in
`metadata$cells` (exact counts) and can be recovered as a raster with
`patch_domain()`. This is deliberate: baking a `-1` sentinel into the label raster
would break every downstream `terra` operation users expect to work.

## 4. Decisions on the explicit questions

**1. Are labeling and metrics separate functions?** Yes.
`label_patches()` is the expensive, cacheable step. `patch_metrics()` is cheap and
users frequently want it several ways (different core depths, different units)
without relabeling. `analyze_patches()` is the convenience wrapper that does both
and is what most workflows should call. Metrics are *not* computed inside
`label_patches()` — but the cell counts are, because they fall out of labeling for
free and are needed for validation.

**2. Does the result hold the raster in memory or on disk?** Both, chosen by the
user. `patch_result$patches` is an active accessor: if the run wrote to disk it
returns `terra::rast(path)`; if it stayed in memory it returns the `SpatRaster` it
already holds. Callers cannot tell the difference, and nothing large is duplicated.
`output = <path>` selects disk; the default is memory.

**3. Background / outside / missing representation?** See §3.

**4. How is perimeter defined?** Perimeter is the summed length of **shared cell
edges** between a patch cell and any cell that is not part of that patch, plus
edges on the grid border. Corner (diagonal) adjacency contributes **zero** — two
cells touching only at a corner share no physical boundary, so counting it would
be geometrically meaningless regardless of the connectivity rule used for
labeling. This means an 8-connected patch's perimeter includes the full outline of
its diagonal "pinch points", which is correct: that is where the physical boundary
runs.

Edge lengths respect non-square cells. An edge shared with a left/right neighbour
runs north–south and has length `yres`; an edge shared with an up/down neighbour
runs east–west and has length `xres`. Getting this backwards is a classic bug, so
the two are accumulated in separate columns (`edge_ns_cells`, `edge_ew_cells`)
that are exact integers and independent of CRS.

**5. Domain edge vs habitat edge?** Split into two reported components that sum to
`perimeter`:

- `edge_valid` — boundary against an in-domain, non-missing cell. This is real
  habitat edge: the patch genuinely stops here.
- `edge_domain` — boundary against an outside-domain cell, a missing cell, or the
  grid border. This is an artefact of the study area, not ecology. The patch may
  well continue; we just cannot see it.

`touches_domain_edge` is `edge_domain > 0`. Any patch with this flag has area and
edge metrics that are **lower bounds**. This is the flag that lets a user drop
boundary-truncated patches from an analysis rather than silently treating a
clipped patch as a small one.

**6. How is a mother patch chosen after a split?** The primary predecessor of a
time-2 patch is the time-1 patch contributing the **largest overlap cell count**.
Ties are broken deterministically: larger time-1 patch first, then lower time-1
ID. Any tie that had to be broken is flagged in `primary_tie`, so it is auditable
rather than invisible.

This rule is a default, not a commitment. The full `overlaps` table is returned
with both directional proportions and both ranks, so a user can re-derive lineage
under any rule they prefer (largest proportion retained, largest area, nearest
centroid) without re-running the expensive step.

**7. How are mergers represented?** In the `overlaps` table as many time-1 rows
mapping to one time-2 patch. The `events` table then summarises per patch with
`n_predecessors` / `n_descendants`, and event classes are assigned from the
connected components of the bipartite overlap graph, so a component with 2 time-1
and 1 time-2 patch is a `merger` and a component with 2-and-2 is `complex`. There
is no attempt to force a merged patch into a single lineage: the component is the
honest unit.

**8. Thresholds on trivial overlaps?** Defaults keep everything
(`min_overlap_cells = 1`, `min_overlap_prop = 0`) because suppression is a
scientific choice and a silent default would be a trap. But event classification
is *acutely* sensitive to slivers — a patch donating one cell to a neighbour would
otherwise register as a split — so both thresholds are exposed and documented, and
`compare_patches()` reports how many links the thresholds dropped.

**9. Geographic (lon/lat) rasters?** Supported, exactly, without pretending
planar geometry. For a lon/lat grid, cell area and east–west extent vary with
latitude but are **constant within a raster row** (verified against
`terra::cellSize()`). So diamondback computes per-row geometry vectors in R —
cell area by row, north–south edge length by row, and east–west edge length at
each of the `nrow+1` horizontal grid lines — and passes them to Python as
`bincount` *weights*. Length and area are then accumulated directly rather than
counted and multiplied afterwards. The same code path serves projected rasters
with constant vectors, so there is no separate lon/lat branch to rot.

**The supported geometry is narrow, and enforced rather than assumed.** The
method is correct for a **north-up, axis-aligned, regular** lon/lat grid, which
is exactly what terra's data model represents: an extent plus a cell count, with
no rotation term. Two things are checked in `db_check_grid()` rather than
trusted:

- **Rotation.** terra can carry a rotated GDAL geotransform, and `is.rotated()`
  reports it. Every cell-index-to-map conversion here (bounding boxes,
  centroids, per-row latitudes) assumes no rotation, so a rotated raster is
  rejected with a pointer to `terra::rectify()` rather than silently producing
  confident nonsense.
- **Latitude bounds.** terra will happily construct a lon/lat raster reaching to
  100°N. The geodesic row geometry is undefined there, so it is an error.

Cells are treated as spherical quadrilaterals bounded by meridians and
parallels. Areas come from `terra::cellSize()`, which is geodesic on the
ellipsoid; edge lengths are geodesic distances at the relevant latitude.

The one place this fails is `patch_core_area()`: the Euclidean distance transform
needs uniform sampling, which a lon/lat grid does not have. That function errors
on lon/lat unless `units = "cells"` is given explicitly.

**10. Memory ceiling and failure behaviour?** Estimated before allocation, from
cell count and dtype, and compared against detected available RAM. Above
`max_memory_frac` (default 0.6) the operation **errors before allocating**, naming
the estimate, the ceiling and the overrides. Peak cost is roughly:

| step | bytes/cell |
|---|---|
| labeling | ~6 (code 1 + bool 1 + int32 labels 4) |
| metrics | ~5 (chunked; no full-array temporaries) |
| core area | ~10 (float64 EDT distances dominate) |

**11. R vs Python?** See §2. One implementation detail is load-bearing enough to
state here: the backend module is imported with **`convert = FALSE`**. With
reticulate's default, every return value is converted to an R object, which
would pull each label array back into R and defeat the entire memory design. Big
arrays therefore stay as Python handles and cross back only through `py_get()`;
small per-patch results are converted explicitly with `py_r()`. This was caught
by a test asserting the arrays were still alive after labelling — the code
"worked" with the default, silently and expensively.

**12. Python exception translation.** Every call crosses `py_try()`, which catches
the condition, extracts the Python exception type and message, and re-raises an R
error of class `diamondback_python_error` naming the operation that failed.
`MemoryError` is special-cased into a message about the memory ceiling rather than
a bare traceback.

**13. What is needed to reproduce a run?** Everything in `metadata`: source
fingerprint (path, size, mtime, and content hash — see §5 for how strong that
is), full geometry (CRS WKT, extent, resolution, dimensions), the mask
fingerprint, `class`, `directions`, `na` handling, cell-state counts, package
version, algorithm version, and Python/NumPy/SciPy versions. This is also exactly
the cache key.

## 4b. Further decisions

**Patch IDs are labels, not identities.** IDs are assigned in raster scan order,
so they are deterministic: identical inputs always give identical numbering. They
are *not* stable across anything else. Cropping, changing the mask, switching
connectivity, or one cell flipping near the top-left renumbers everything after
it, and patch 47 in 1985 has nothing to do with patch 47 in 2024. Nothing in the
API invites a join on `patch_id` across runs, and the documentation says so
outright, because it is the obvious thing for a user to try and it is silently
wrong. Correspondence between runs comes only from `compare_patches()`, derived
from actual cell overlap. `lineage_id` is likewise scoped to one comparison.

**Threading and determinism.** Version 1 is single-process and single-threaded,
and the results are bit-for-bit reproducible. `scipy.ndimage.label()`,
`distance_transform_edt()`, `find_objects()` and the NumPy reductions used here
(`bincount`, `unique`) are all serial C. None of them dispatch to BLAS, so
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and friends have no effect on either
the result or the runtime — there is no thread-count-dependent reduction order to
perturb the arithmetic. Parallelism would have to come from a different labelling
backend (§6), not from tuning environment variables.

**The class limit.** Classes are encoded as `3 + k` in the cell-state array, so
the ceiling is the dtype: 253 classes in `uint8`, 65,533 in `uint16`. Version 1
picks the narrower dtype when it fits and widens automatically when it does not,
because 252 was an arbitrary artefact of an encoding choice and real categorical
rasters exceed it. The dtype is chosen once per run from the class count — not
per block from the values a block happens to contain — and every comparison is by
value, so nothing downstream knows or cares which is in use. `uint16` costs one
extra byte per cell, which the memory estimate accounts for.

Multi-class runs label each class in a **separate pass** over the raster, so cost
scales with the class count; past 100 classes `class = "all"` warns rather than
quietly taking an hour.

## 5. Caching

Cache identity is a SHA-1 of the *semantic* inputs: source fingerprint, geometry,
class, mask fingerprint, connectivity, NA policy, and `algorithm_version`. Note
what is deliberately **excluded**: `output` paths, progress settings, and metric
options do not invalidate a labeling cache, because they cannot change the labels.

`algorithm_version` is a package constant bumped by hand whenever the labeling
kernel changes semantics. Package version alone is too coarse (a typo fix in a
docstring should not throw away an hour of labeling) and too fine at once.

Stale results are never reused silently: a fingerprint mismatch is reported with
the specific field that differs, then the work is redone.

**How strong is the fingerprint?** `fingerprint = "auto"` (the default) hashes
file contents under 200 MB and falls back to path + size + mtime above that,
because hashing a 40 GB raster to decide whether to reuse a cache can cost more
than the labelling did. Size and mtime are a *heuristic*, not proof: a file
rewritten at the same size with its timestamp restored would be treated as
unchanged. `fingerprint = "full"` therefore always hashes, however large, and
`"fast"` never does. The mode is part of the cache key, so a result cached under
a weak fingerprint is never silently reused for a request that asked for a strong
one — and the mode is the *first* field compared, so switching modes reports
"fingerprint differs" rather than the technically-true-but-misleading "hash
differs".

**Canonical fields, not live values.** One side of a cache comparison has been
through JSON and the other has not, and JSON does not preserve R's type
distinctions: `terra::ext()` returns *named* numerics, `0.0` comes back as an
integer, and `NA_real_` does not come back as `NA_real_`. Comparing the live
values produced a cache that never hit — which is worse than it sounds, because
it fails silently and just costs an hour per run. So every cache field is
canonicalised to a character scalar (`db_cache_fields()`), and those strings are
written into `metadata.json` and compared on the way back in. Numbers are
formatted at 15 significant digits and the JSON is written at full double
precision, so a lon/lat resolution of 1/12° survives the round trip exactly.

An in-memory raster larger than a million cells has no file to fingerprint and
no cheap content hash, so it is **never** served from cache. A slow rerun beats
a wrong answer.

## 6. Large rasters

Version 1 requires the label array to fit in RAM. This is not laziness, it is
`scipy.ndimage.label()`'s contract.

What is done anyway to make that ceiling as high as possible:

- **Blockwise ingestion.** The source is read in row blocks and written straight
  into a preallocated `uint8` NumPy array. The whole raster is never an R numeric
  matrix. This alone is an 8× reduction against the obvious `values(x)` approach.
- **Tight cropping.** With `crop = TRUE` (default when a mask is given) the grid
  is cropped to the mask's extent before allocation.
- **Chunked reductions.** Metrics, edges and overlaps run in row chunks, so no
  step allocates a full-array float temporary.
- **dtype selection.** Labels are `int32` unless the cell count could overflow it,
  in which case `int64`.

What was considered and rejected for v1: Dask and tiled union-find labeling. Both
are real answers to rasters that genuinely exceed RAM, and both are a large amount
of machinery to get right (cross-tile label merging is where these go wrong). The
in-memory path is correct and fast for the target case. The labeling kernel is
isolated behind `label_array()` so a second backend can be added without touching
the R API.

Temp files: diamondback tracks every file it creates, **with its provenance**,
and removes only its own temporaries. Ownership is recorded at creation, not
inferred from the path: a user's `output` can perfectly well live under
`tempdir()`, and that makes it a temporary *location*, not a file we are
entitled to delete. `terra::tmpFiles(remove = TRUE)` is never called — it would
invalidate live `SpatRaster` objects belonging to the user, which is precisely
the failure this package exists to end.

## 7. Python environment

**Running an analysis never downloads or installs anything.** This is a hard
rule, and it took a correction to get right.

The obvious design is `reticulate::py_require()` (reticulate >= 1.41): declare
numpy and scipy, and let reticulate resolve a uv-managed environment on first
use. It is genuinely convenient and it does not touch the user's existing
environments. But on a machine without a suitable Python it will fetch a whole
CPython interpreter plus numpy and scipy — ~50 MB of downloads, ~336 MB on disk
— triggered by nothing more than a call to `label_patches()`. "It doesn't modify
your environments, it creates a new one" is a dodge. A call to an analysis
function is not consent to download 300 MB.

So `py_require()` is only ever called after explicit consent.
`db_python_mode()` decides where Python comes from **without initialising an
interpreter**, and resolves in this order:

| order | source | installs? |
|---|---|---|
| 1 | Python already initialised in the session | no |
| 2 | `RETICULATE_PYTHON`, or `options(diamondback.python = <path>)` | no |
| 3 | the managed environment, **if** `diamondback_install_python()` was run | already fetched |
| 4 | any Python on `PATH` / in an active venv or conda env that already has numpy+scipy | no |
| 5 | nothing found → **error with instructions** | no |

Step 4 probes by shelling out (`python -c "import numpy, scipy"`) rather than
going through reticulate, because initialising reticulate is itself what can
trigger provisioning. Step 5 is the whole point: diamondback stops and tells the
user their options rather than deciding for them.

`diamondback_install_python()` is the only thing in the package that downloads.
It states what it will fetch and roughly how large it is, asks for confirmation
when interactive, and records consent in
`tools::R_user_dir("diamondback", "config")` so later sessions do not re-ask.
`diamondback_remove_python()` revokes it. The package's own test suite opts in
explicitly (`options(diamondback.python = "managed")` in a test helper) — that
is a deliberate, visible choice in one place, not a default.

`diamondback_check()` is the single diagnostic and reports which of the five
sources was used; `diamondback_python()` is the single place any Python is
touched.

## 8. Validation

`validate_patch_result()` is the collected correctness checks: label consecutiveness,
foreground count vs summed patch sizes, table vs raster patch counts, no foreground
outside the mask, ID range vs integer limits, and metric-internal consistency
(edge components summing to perimeter, core area not exceeding total area). It runs
automatically when `validate = TRUE` (the default is `FALSE` for large runs; the
cost is a second pass).

The test suite independently checks SciPy results against `terra::patches()` on
small rasters under both connectivity rules, at foreground densities from 5% to
95%, and with `NA` present. The two must agree exactly on patch count *and* on
the induced partition of cells, not merely the count — and they do.

**Testing the chunked paths.** The reductions run in 4-million-cell row chunks,
which no test raster would ever reach, so the seam-crossing logic would be
exercised only in production. `CHUNK_CELLS` is therefore a module global rather
than a default argument, and the tests shrink it to force real seams through
30-row rasters. A patch spanning a seam must measure identically to one that
does not. The lon/lat case matters most here: per-row edge lengths are indexed
by absolute row, so an off-by-one at a chunk boundary would silently shift a
patch's latitude.
