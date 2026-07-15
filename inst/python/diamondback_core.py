"""Array kernels for the diamondback R package.

This module is loaded by R via reticulate.import_from_path(); it is not meant to
be run standalone and deliberately has no CLI.

Cell state codes (uint8), see DESIGN.md section 3:
    0      outside the analysis domain
    1      missing / unknown (NA in source, inside the mask)
    2      valid background
    3 + k  foreground of class k

Conventions
-----------
* Arrays are 2D, C-order, shape (nrow, ncol); row 0 is the top row of the raster.
* "vertical edge"   = boundary shared with a left/right neighbour; it runs
                      north-south and its length is the cell height.
* "horizontal edge" = boundary shared with an up/down neighbour; it runs
                      east-west and its length is the cell width.
* Per-label result arrays are length n + 1 and indexed by label; element 0 is
  background and is not meaningful.
"""

import numpy as np
from scipy import ndimage as ndi

__all__ = [
    "versions",
    "code_block",
    "code_counts",
    "label_array",
    "combine_labels",
    "set_rows",
    "patch_stats",
    "patch_bboxes",
    "edge_lengths",
    "core_counts",
    "overlap_counts",
    "output_rows",
    "domain_rows",
    "foreground_outside_domain",
]

ALGORITHM_VERSION = "1"


def versions():
    """Report backend versions; used by diamondback_check()."""
    import scipy
    import sys

    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "algorithm_version": ALGORITHM_VERSION,
    }


# Target cells per row-chunk in the reduction passes. A module global rather
# than a default argument so tests can shrink it and force real chunk seams
# through small rasters -- a patch spanning a seam must measure identically to
# one that does not, and that is not something to leave untested.
CHUNK_CELLS = 4_000_000


def _chunk_bounds(nrow, ncol, target_cells=None):
    """Row-chunk size keeping working arrays bounded regardless of raster width."""
    if target_cells is None:
        target_cells = CHUNK_CELLS
    step = max(1, int(target_cells // max(1, ncol)))
    return [(r, min(r + step, nrow)) for r in range(0, nrow, step)]


def _structure(directions):
    if int(directions) == 4:
        return ndi.generate_binary_structure(2, 1)
    if int(directions) == 8:
        return ndi.generate_binary_structure(2, 2)
    raise ValueError("directions must be 4 or 8, got %r" % (directions,))


# --------------------------------------------------------------------------
# Building the code array
# --------------------------------------------------------------------------

def code_block(vals, nrow, ncol, mask=None, class_values=None):
    """Turn one block of raster values into cell state codes.

    Parameters
    ----------
    vals : 1D array-like, row-major, length nrow*ncol. NaN marks NA.
    mask : 1D array-like or None. Non-zero and non-NaN means "inside domain".
    class_values : sequence of numbers, or None for binary treatment where any
        non-zero value is foreground.

    Precedence is outside-domain > missing > class, so an NA cell outside the
    mask is coded 0.
    """
    v = np.asarray(vals, dtype=np.float64).reshape(int(nrow), int(ncol))
    out = np.full(v.shape, 2, dtype=np.uint8)

    if class_values is None:
        out[v != 0] = 3
    else:
        cv = np.atleast_1d(np.asarray(class_values, dtype=np.float64))
        if cv.size > 252:
            raise ValueError("at most 252 classes are supported in one run")
        for k in range(cv.size):
            out[v == cv[k]] = 3 + k

    out[np.isnan(v)] = 1

    if mask is not None:
        m = np.asarray(mask, dtype=np.float64).reshape(int(nrow), int(ncol))
        out[np.isnan(m) | (m == 0)] = 0

    return out


def code_counts(code):
    """Count cells in each state. Returns a length-256 int64 array."""
    return np.bincount(np.asarray(code).ravel(), minlength=256).astype(np.int64)


def foreground_outside_domain(code, class_index):
    """Always zero by construction; recomputed as an independent check."""
    fg = np.asarray(code) == (3 + int(class_index))
    return int(np.count_nonzero(fg & (np.asarray(code) == 0)))


# --------------------------------------------------------------------------
# Labeling
# --------------------------------------------------------------------------

def label_array(code, class_index, directions=8, use_int64=False):
    """Connected-component labeling of one class.

    Returns a dict {"labels", "n"} rather than a tuple: the R side holds the
    module with convert=False so that big arrays stay in Python, and named dict
    access is far less error-prone than 0-based tuple indexing across the
    boundary.

    Labels are 0 for every non-foreground cell and 1..n within the class.
    """
    code = np.asarray(code)
    fg = code == (3 + int(class_index))
    dtype = np.int64 if use_int64 else np.int32
    labels, n = ndi.label(fg, structure=_structure(directions), output=dtype)
    return {"labels": labels, "n": int(n)}


def set_rows(arr, r0, r1, block):
    """Write a row block into a preallocated array, in place.

    Doing this in Python rather than through reticulate's slice assignment
    keeps the ingestion path explicit and avoids materialising the block on the
    R side a second time.
    """
    arr[int(r0):int(r1)] = block
    return None


def combine_labels(labels_list, offsets, shape, use_int64=False):
    """Merge per-class label arrays into one array with globally unique IDs."""
    dtype = np.int64 if use_int64 else np.int32
    out = np.zeros(tuple(shape), dtype=dtype)
    for lab, off in zip(labels_list, offsets):
        nz = lab > 0
        out[nz] = lab[nz] + int(off)
    return out


# --------------------------------------------------------------------------
# Per-patch reductions
# --------------------------------------------------------------------------

def patch_stats(labels, n, cell_area_by_row=None):
    """Cell count, centroid accumulators and (optionally) area, in one pass.

    cell_area_by_row : length-nrow array of per-cell area, used for lon/lat
        grids where area varies with latitude. None means "count only".

    Returns dict of length n+1 arrays: count, sum_row, sum_col, area.
    Centroid sums are in 0-based cell index space (row/col centre indices).
    """
    labels = np.asarray(labels)
    n = int(n)
    nrow, ncol = labels.shape

    count = np.zeros(n + 1, dtype=np.int64)
    sum_row = np.zeros(n + 1, dtype=np.float64)
    sum_col = np.zeros(n + 1, dtype=np.float64)
    area = np.zeros(n + 1, dtype=np.float64) if cell_area_by_row is not None else None
    if cell_area_by_row is not None:
        cell_area_by_row = np.asarray(cell_area_by_row, dtype=np.float64)

    cols = np.arange(ncol, dtype=np.float64)

    for r0, r1 in _chunk_bounds(nrow, ncol):
        lab = labels[r0:r1].ravel()
        sel = lab > 0
        if not sel.any():
            continue
        l = lab[sel]
        count += np.bincount(l, minlength=n + 1)

        rr = np.repeat(np.arange(r0, r1, dtype=np.float64), ncol)[sel]
        cc = np.tile(cols, r1 - r0)[sel]
        sum_row += np.bincount(l, weights=rr, minlength=n + 1)
        sum_col += np.bincount(l, weights=cc, minlength=n + 1)

        if area is not None:
            aa = np.repeat(cell_area_by_row[r0:r1], ncol)[sel]
            area += np.bincount(l, weights=aa, minlength=n + 1)

    return {"count": count, "sum_row": sum_row, "sum_col": sum_col, "area": area}


def patch_bboxes(labels, n):
    """Bounding boxes as 0-based inclusive cell indices.

    Returns dict of length n+1 arrays: row_min, row_max, col_min, col_max.
    Uses scipy's find_objects, which is a single C pass over the array.
    """
    n = int(n)
    objs = ndi.find_objects(np.asarray(labels), max_label=n)
    row_min = np.full(n + 1, -1, dtype=np.int64)
    row_max = np.full(n + 1, -1, dtype=np.int64)
    col_min = np.full(n + 1, -1, dtype=np.int64)
    col_max = np.full(n + 1, -1, dtype=np.int64)
    for i, sl in enumerate(objs):
        if sl is None:
            continue
        lab = i + 1
        row_min[lab] = sl[0].start
        row_max[lab] = sl[0].stop - 1
        col_min[lab] = sl[1].start
        col_max[lab] = sl[1].stop - 1
    return {
        "row_min": row_min,
        "row_max": row_max,
        "col_min": col_min,
        "col_max": col_max,
    }


# --------------------------------------------------------------------------
# Boundary edges
# --------------------------------------------------------------------------

def edge_lengths(labels, code, n, dy_by_row, dx_by_gridline, na_background=False):
    """Accumulate exposed patch boundary, split by what is on the other side.

    An edge is counted when a patch cell is adjacent (rook only) to a cell with
    a different label, or lies on the grid border. Diagonal contact contributes
    nothing: corner-touching cells share no physical boundary.

    "valid"  = the neighbour is in-domain and not missing (real habitat edge).
    "domain" = the neighbour is outside the domain, missing, or off-grid
               (an artefact of the study area).

    dy_by_row : length-nrow, north-south extent of a cell in each row. Weighs
        vertical edges (left/right neighbours).
    dx_by_gridline : length-nrow+1, east-west extent at each horizontal grid
        line. Weighs horizontal edges (up/down neighbours); element j is the
        line above row j.

    Returns dict of length n+1 arrays. The *_cells entries are exact integer
    edge counts and carry no CRS assumption; the *_len entries are map units.
    """
    labels = np.asarray(labels)
    code = np.asarray(code)
    n = int(n)
    nrow, ncol = labels.shape
    dy = np.asarray(dy_by_row, dtype=np.float64)
    dx = np.asarray(dx_by_gridline, dtype=np.float64)

    valid_len = np.zeros(n + 1, dtype=np.float64)
    domain_len = np.zeros(n + 1, dtype=np.float64)
    ns_cells = np.zeros(n + 1, dtype=np.int64)   # vertical edges
    ew_cells = np.zeros(n + 1, dtype=np.int64)   # horizontal edges
    valid_cells = np.zeros(n + 1, dtype=np.int64)
    domain_cells = np.zeros(n + 1, dtype=np.int64)

    def _is_valid(c):
        # A neighbour counts as real habitat edge if it is in-domain and known.
        if na_background:
            return c >= 1
        return c >= 2

    def _accum(lab_sel, nb_valid, w, vertical):
        """lab_sel: labels of exposed cells; nb_valid: bool, same shape."""
        if lab_sel.size == 0:
            return
        vs = nb_valid
        if vs.any():
            l = lab_sel[vs]
            valid_len[:] += np.bincount(l, weights=w[vs], minlength=n + 1)
            valid_cells[:] += np.bincount(l, minlength=n + 1)
        ds = ~nb_valid
        if ds.any():
            l = lab_sel[ds]
            domain_len[:] += np.bincount(l, weights=w[ds], minlength=n + 1)
            domain_cells[:] += np.bincount(l, minlength=n + 1)
        tgt = ns_cells if vertical else ew_cells
        tgt[:] += np.bincount(lab_sel, minlength=n + 1)

    for r0, r1 in _chunk_bounds(nrow, ncol):
        # ---- vertical edges: left/right neighbours within each row ----
        L = labels[r0:r1, :-1]
        R = labels[r0:r1, 1:]
        if L.size:
            diff = L != R
            wy = np.broadcast_to(dy[r0:r1, None], L.shape)
            vR = _is_valid(code[r0:r1, 1:])
            s = diff & (L > 0)
            _accum(L[s], vR[s], wy[s], True)
            vL = _is_valid(code[r0:r1, :-1])
            s = diff & (R > 0)
            _accum(R[s], vL[s], wy[s], True)

        # ---- vertical edges on the left and right grid borders ----
        # For a single-column raster this visits column 0 twice, which is
        # correct: such a cell is exposed on both its left and right side.
        for col in (0, ncol - 1):
            lab = labels[r0:r1, col]
            s = lab > 0
            if s.any():
                _accum(lab[s], np.zeros(int(s.sum()), dtype=bool), dy[r0:r1][s], True)

        # ---- horizontal edges: up/down neighbours, pairs (j, j+1) ----
        j1 = min(r1, nrow - 1)
        if j1 > r0:
            U = labels[r0:j1]
            D = labels[r0 + 1:j1 + 1]
            diff = U != D
            wx = np.broadcast_to(dx[r0 + 1:j1 + 1, None], U.shape)
            vD = _is_valid(code[r0 + 1:j1 + 1])
            s = diff & (U > 0)
            _accum(U[s], vD[s], wx[s], False)
            vU = _is_valid(code[r0:j1])
            s = diff & (D > 0)
            _accum(D[s], vU[s], wx[s], False)

    # ---- horizontal edges on the top and bottom grid borders ----
    # As above, a single-row raster visits row 0 twice, once per border.
    for row, gl in ((0, 0), (nrow - 1, nrow)):
        lab = labels[row]
        s = lab > 0
        if s.any():
            k = int(s.sum())
            _accum(lab[s], np.zeros(k, dtype=bool), np.full(k, dx[gl]), False)

    return {
        "edge_valid_len": valid_len,
        "edge_domain_len": domain_len,
        "edge_ns_cells": ns_cells,
        "edge_ew_cells": ew_cells,
        "edge_valid_cells": valid_cells,
        "edge_domain_cells": domain_cells,
    }


# --------------------------------------------------------------------------
# Core area
# --------------------------------------------------------------------------

def core_counts(code, labels, n, class_index, depth, sampling, edge_mode="all"):
    """Cells further than `depth` from an edge, per patch, via an exact EDT.

    Distance is centre-to-centre: the transform returns, for each foreground
    cell, the distance to the centre of the nearest edge-source cell. A cell
    immediately adjacent to non-habitat therefore has distance = one cell width,
    not zero. Core cells are those with distance strictly greater than `depth`.

    edge_mode "all"        : every non-foreground cell is an edge source, and
                             the area beyond the grid is treated as non-habitat.
    edge_mode "background" : only valid background is an edge source. Missing
                             and outside-domain cells, and the area beyond the
                             grid, are treated as habitat, so patches clipped by
                             the study boundary are not penalised.
    """
    code = np.asarray(code)
    labels = np.asarray(labels)
    n = int(n)
    fg = code == (3 + int(class_index))

    if edge_mode == "all":
        src = fg
        pad_value = False          # off-grid is non-habitat -> an edge
    elif edge_mode == "background":
        src = code != 2            # zeros are valid background only
        pad_value = True           # off-grid is not an edge
    else:
        raise ValueError("edge_mode must be 'all' or 'background', got %r" % (edge_mode,))

    padded = np.pad(src, 1, mode="constant", constant_values=pad_value)
    dist = ndi.distance_transform_edt(padded, sampling=tuple(float(s) for s in sampling))
    dist = dist[1:-1, 1:-1]

    core = fg & (dist > float(depth))
    del dist, padded, src

    counts = np.bincount(labels[core], minlength=n + 1).astype(np.int64)
    counts[0] = 0
    return {"core_count": counts, "core_mask": core}


# --------------------------------------------------------------------------
# Temporal overlap
# --------------------------------------------------------------------------

def overlap_counts(lab1, lab2, n1, n2):
    """Cross-tabulate overlapping cells between two label arrays.

    Only cells that are inside a patch in both arrays are counted. Returns
    parallel arrays (id1, id2, cells). Chunked so no full-array key vector is
    materialised.
    """
    lab1 = np.asarray(lab1)
    lab2 = np.asarray(lab2)
    if lab1.shape != lab2.shape:
        raise ValueError("label arrays have different shapes: %r vs %r"
                         % (lab1.shape, lab2.shape))
    n1, n2 = int(n1), int(n2)

    stride = n2 + 1
    if (n1 + 1) * stride > 2 ** 62:
        raise OverflowError(
            "patch counts too large to cross-tabulate: %d x %d" % (n1, n2))

    nrow, ncol = lab1.shape
    keys_parts = []
    counts_parts = []
    for r0, r1 in _chunk_bounds(nrow, ncol):
        a = lab1[r0:r1].ravel()
        b = lab2[r0:r1].ravel()
        sel = (a > 0) & (b > 0)
        if not sel.any():
            continue
        k = a[sel].astype(np.int64) * stride + b[sel].astype(np.int64)
        uk, uc = np.unique(k, return_counts=True)
        keys_parts.append(uk)
        counts_parts.append(uc)

    if not keys_parts:
        return {
            "id1": np.zeros(0, dtype=np.int64),
            "id2": np.zeros(0, dtype=np.int64),
            "cells": np.zeros(0, dtype=np.int64),
        }

    keys = np.concatenate(keys_parts)
    counts = np.concatenate(counts_parts).astype(np.int64)
    uk, inv = np.unique(keys, return_inverse=True)
    tot = np.bincount(inv, weights=counts, minlength=uk.size).astype(np.int64)

    return {"id1": uk // stride, "id2": uk % stride, "cells": tot}


# --------------------------------------------------------------------------
# Handing labels back to R
# --------------------------------------------------------------------------

INT32_NA = np.int32(-2147483648)  # R's NA_integer_


def output_rows(labels, code, r0, r1, na_background=False):
    """One row block of the labelled raster, as int32 with R's NA sentinel.

    Cells outside the domain (and missing cells, unless na_background) become
    NA. reticulate maps int32 INT_MIN straight onto R's NA_integer_, so no
    second pass is needed on the R side.
    """
    lab = np.asarray(labels)[int(r0):int(r1)]
    cod = np.asarray(code)[int(r0):int(r1)]
    out = lab.astype(np.int32, copy=True)
    if na_background:
        bad = cod == 0
    else:
        bad = cod <= 1
    out[bad] = INT32_NA
    return out.ravel()


def domain_rows(code, r0, r1):
    """One row block of the cell-state code array, for patch_domain()."""
    return np.asarray(code)[int(r0):int(r1)].ravel().astype(np.int32)
