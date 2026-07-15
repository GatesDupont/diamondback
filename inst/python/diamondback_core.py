"""Array kernels for the diamondback R package.

This module is loaded by R via reticulate.import_from_path(); it is not meant to
be run standalone and deliberately has no CLI.

Cell state codes, see DESIGN.md section 3:
    0      outside the analysis domain
    1      missing / unknown (NA in source, inside the mask)
    2      valid background
    3 + k  foreground of class k

The code array is uint8 (1 byte/cell) when there are few enough classes to fit,
and uint16 (2 bytes/cell) otherwise. The dtype is an encoding detail: every
comparison is by value, so nothing downstream needs to know which is in use.

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

import math

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

# Labeling semantics. Changing this invalidates caches. See DESIGN.md section 5.
ALGORITHM_VERSION = "2"

# The R-facing surface of this module: function names, their signatures, and the
# keys of the dicts they return. Bump it whenever any of those change, and bump
# PY_INTERFACE_VERSION in R/python.R to match.
#
# This exists because the R code and this file are two halves of one program that
# can drift apart -- most easily when a package is reinstalled while an R session
# still holds the old namespace in memory. The R half then calls the new Python
# half with the old expectations and gets, say, KeyError: 'edge_domain_len',
# which tells the user nothing. The handshake turns that into a sentence.
INTERFACE_VERSION = "3"


def versions():
    """Report backend versions; used by diamondback_check() and the handshake."""
    import scipy
    import sys

    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "algorithm_version": ALGORITHM_VERSION,
        "interface_version": INTERFACE_VERSION,
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

# Codes 0-2 are the non-class states, so class k occupies code 3 + k.
MAX_CLASSES_UINT8 = 253    # codes 3..255
MAX_CLASSES_UINT16 = 65533  # codes 3..65535


def code_dtype(n_classes):
    """Smallest unsigned dtype that can encode this many classes.

    Returns an np.dtype instance rather than the scalar type, so callers can
    inspect .name and pass it straight to np.zeros().
    """
    n = int(n_classes)
    if n <= MAX_CLASSES_UINT8:
        return np.dtype(np.uint8)
    if n <= MAX_CLASSES_UINT16:
        return np.dtype(np.uint16)
    raise ValueError(
        "at most %d classes are supported in one run, got %d"
        % (MAX_CLASSES_UINT16, n))


def code_block(vals, nrow, ncol, mask=None, class_groups=None, n_classes=1):
    """Turn one block of raster values into cell state codes.

    Parameters
    ----------
    vals : 1D array-like, row-major, length nrow*ncol. NaN marks NA.
    mask : 1D array-like or None. Non-zero and non-NaN means "inside domain".
    class_groups : sequence of groups, each group a sequence of raster values,
        or None for binary treatment where any non-zero value is foreground.
        Every value in group k is coded 3 + k, so values sharing a group become
        one class and are labelled as continuous habitat across their
        boundaries; values in different groups get different codes and are
        labelled separately.
    n_classes : number of groups in the run, which fixes the dtype. Passed
        explicitly so that every block of a run agrees, rather than each block
        choosing from the values it happens to contain.

    Precedence is outside-domain > missing > class, so an NA cell outside the
    mask is coded 0.
    """
    v = np.asarray(vals, dtype=np.float64).reshape(int(nrow), int(ncol))
    dt = code_dtype(n_classes)
    out = np.full(v.shape, 2, dtype=dt)

    if class_groups is None:
        out[v != 0] = 3
    else:
        code_dtype(len(class_groups))  # validates the count
        for k, grp in enumerate(class_groups):
            g = np.atleast_1d(np.asarray(grp, dtype=np.float64))
            # One isin() per group: every value in the group lands on the same
            # code, which is what makes a grouped class one connected surface.
            out[np.isin(v, g)] = 3 + k

    out[np.isnan(v)] = 1

    if mask is not None:
        m = np.asarray(mask, dtype=np.float64).reshape(int(nrow), int(ncol))
        out[np.isnan(m) | (m == 0)] = 0

    return out


def code_counts(code):
    """Count cells in each state, indexed by code.

    Length follows the dtype, so callers must not assume 256.
    """
    code = np.asarray(code)
    n = 256 if code.dtype == np.uint8 else int(MAX_CLASSES_UINT16) + 3
    return np.bincount(code.ravel(), minlength=n).astype(np.int64)


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

    The far side is classified into three kinds, mirroring the cell states.
    Lumping the last two together would say that a patch beside an unsurveyed
    hole and a patch cut off by the study-area boundary are the same thing, and
    they are not:

    "valid"   = an in-domain, known cell (background, or another class).
                Real habitat edge: the patch genuinely stops here.
    "missing" = a cell inside the domain whose value is unknown (code 1).
                The patch may continue into it; nobody looked.
    "outside" = a cell excluded by the mask (code 0), or off the grid.
                An artefact of where the study area was drawn.

    With na_background, the caller has declared NA to mean genuine absence, so
    missing cells count as real habitat edge and "missing" is always zero.

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
    missing_len = np.zeros(n + 1, dtype=np.float64)
    outside_len = np.zeros(n + 1, dtype=np.float64)
    ns_cells = np.zeros(n + 1, dtype=np.int64)   # vertical edges
    ew_cells = np.zeros(n + 1, dtype=np.int64)   # horizontal edges
    valid_cells = np.zeros(n + 1, dtype=np.int64)
    missing_cells = np.zeros(n + 1, dtype=np.int64)
    outside_cells = np.zeros(n + 1, dtype=np.int64)

    # Neighbour kind: 0 = outside, 1 = missing, 2 = valid. With na_background a
    # missing neighbour is real habitat edge, so it is promoted to valid.
    def _nb_kind(c):
        k = np.where(c >= 2, 2, c).astype(np.uint8)
        if na_background:
            k = np.where(k == 1, 2, k).astype(np.uint8)
        return k

    def _accum(lab_sel, kind, w, vertical):
        """lab_sel: labels of exposed cells; kind: neighbour kind, same shape."""
        if lab_sel.size == 0:
            return
        for k, (lsum, lcnt) in enumerate((
            (outside_len, outside_cells),
            (missing_len, missing_cells),
            (valid_len, valid_cells),
        )):
            sel = kind == k
            if sel.any():
                l = lab_sel[sel]
                lsum[:] += np.bincount(l, weights=w[sel], minlength=n + 1)
                lcnt[:] += np.bincount(l, minlength=n + 1)
        tgt = ns_cells if vertical else ew_cells
        tgt[:] += np.bincount(lab_sel, minlength=n + 1)

    for r0, r1 in _chunk_bounds(nrow, ncol):
        # ---- vertical edges: left/right neighbours within each row ----
        L = labels[r0:r1, :-1]
        R = labels[r0:r1, 1:]
        if L.size:
            diff = L != R
            wy = np.broadcast_to(dy[r0:r1, None], L.shape)
            kR = _nb_kind(code[r0:r1, 1:])
            s = diff & (L > 0)
            _accum(L[s], kR[s], wy[s], True)
            kL = _nb_kind(code[r0:r1, :-1])
            s = diff & (R > 0)
            _accum(R[s], kL[s], wy[s], True)

        # ---- vertical edges on the left and right grid borders ----
        # For a single-column raster this visits column 0 twice, which is
        # correct: such a cell is exposed on both its left and right side.
        # Off-grid is "outside": the study area stops here.
        for col in (0, ncol - 1):
            lab = labels[r0:r1, col]
            s = lab > 0
            if s.any():
                _accum(lab[s], np.zeros(int(s.sum()), dtype=np.uint8), dy[r0:r1][s], True)

        # ---- horizontal edges: up/down neighbours, pairs (j, j+1) ----
        j1 = min(r1, nrow - 1)
        if j1 > r0:
            U = labels[r0:j1]
            D = labels[r0 + 1:j1 + 1]
            diff = U != D
            wx = np.broadcast_to(dx[r0 + 1:j1 + 1, None], U.shape)
            kD = _nb_kind(code[r0 + 1:j1 + 1])
            s = diff & (U > 0)
            _accum(U[s], kD[s], wx[s], False)
            kU = _nb_kind(code[r0:j1])
            s = diff & (D > 0)
            _accum(D[s], kU[s], wx[s], False)

    # ---- horizontal edges on the top and bottom grid borders ----
    # As above, a single-row raster visits row 0 twice, once per border.
    for row, gl in ((0, 0), (nrow - 1, nrow)):
        lab = labels[row]
        s = lab > 0
        if s.any():
            k = int(s.sum())
            _accum(lab[s], np.zeros(k, dtype=np.uint8), np.full(k, dx[gl]), False)

    return {
        "edge_valid_len": valid_len,
        "edge_missing_len": missing_len,
        "edge_outside_len": outside_len,
        "edge_ns_cells": ns_cells,
        "edge_ew_cells": ew_cells,
        "edge_valid_cells": valid_cells,
        "edge_missing_cells": missing_cells,
        "edge_outside_cells": outside_cells,
    }


# --------------------------------------------------------------------------
# Core area
# --------------------------------------------------------------------------

def _edt_row_strips(nrow, ncol, halo_rows, target_bytes):
    """Row strips for the tiled distance transform, sized to a memory target.

    Strips span the full width, so no column halo is ever needed and both
    vertical edges of every strip are true grid borders. That is a real
    simplification, not a shortcut: it removes an entire class of tile-seam
    bookkeeping.

    Yields (r0, r1) tile bounds; the halo is added by the caller.
    """
    # bool src + float64 dist, over the padded strip
    per_row = ncol * 9
    h = int(target_bytes // per_row) - 2 * halo_rows
    h = max(1, h)
    if h >= nrow:
        return [(0, nrow)]
    return [(r, min(r + h, nrow)) for r in range(0, nrow, h)]


def core_counts(code, labels, n, class_index, depth, sampling, edge_mode="all",
                want_mask=False, target_bytes=256 * 1024 * 1024):
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

    Tiling
    ------
    A full-array transform holds a float64 distance for every cell -- 8 bytes,
    which on a 724M-cell raster is 5.8 GB and simply does not fit on an ordinary
    machine. So the work is done in row strips, and the result is *identical*
    rather than approximate.

    Why it is exact. Give each strip a halo of ceil(depth / cell size) rows.
    Then for any cell in the strip interior:

      * if its true distance to a source is <= depth, that source lies within
        `depth` of it, hence within the halo, hence inside the window -- so the
        windowed transform finds it and returns the true distance;
      * if its true distance is > depth, the windowed transform can only ever
        report a *larger* value (it sees a subset of the sources), so it is
        still > depth.

    Either way the test `dist > depth` gives the same answer as the full-array
    transform. This holds only because the test is a threshold at `depth`; the
    distances themselves are not trustworthy beyond the halo, and are never
    returned.

    The other half of the argument is padding. Only sides that are a true grid
    border are padded, and by exactly the one ring the full-array version used.
    Interior strip edges get real data as their halo and are never padded --
    padding them would invent sources that do not exist and report cells as
    edge-affected when they are not.
    """
    code = np.asarray(code)
    labels = np.asarray(labels)
    n = int(n)
    depth = float(depth)
    nrow, ncol = code.shape
    sy, sx = float(sampling[0]), float(sampling[1])
    fg_code = 3 + int(class_index)

    if edge_mode == "all":
        pad_value = False          # off-grid is non-habitat -> an edge
    elif edge_mode == "background":
        pad_value = True           # off-grid is not an edge
    else:
        raise ValueError("edge_mode must be 'all' or 'background', got %r" % (edge_mode,))

    # A source within Euclidean `depth` is at most depth/sy rows away, so that
    # many rows of halo is sufficient -- and sufficiency is the whole proof.
    halo = int(math.ceil(depth / sy)) if sy > 0 else 0

    counts = np.zeros(n + 1, dtype=np.int64)
    mask_out = np.zeros((nrow, ncol), dtype=bool) if want_mask else None

    for r0, r1 in _edt_row_strips(nrow, ncol, halo, target_bytes):
        w0, w1 = max(0, r0 - halo), min(nrow, r1 + halo)
        sub = code[w0:w1]
        src = (sub == fg_code) if edge_mode == "all" else (sub != 2)

        # Pad only true borders. Strips span the full width, so left and right
        # always are; top and bottom only for the first and last strip.
        pad_top = 1 if w0 == 0 else 0
        pad_bot = 1 if w1 == nrow else 0
        src = np.pad(src, ((pad_top, pad_bot), (1, 1)),
                     mode="constant", constant_values=pad_value)

        dist = ndi.distance_transform_edt(src, sampling=(sy, sx))
        dist = dist[pad_top:dist.shape[0] - pad_bot, 1:-1]

        # Drop the halo, keeping the strip proper.
        d = dist[(r0 - w0):(r0 - w0) + (r1 - r0)]
        fg = code[r0:r1] == fg_code
        core = fg & (d > depth)
        del dist, d, src, fg

        if core.any():
            counts += np.bincount(labels[r0:r1][core], minlength=n + 1)
        if want_mask:
            mask_out[r0:r1] = core

    counts[0] = 0
    return {"core_count": counts.astype(np.int64), "core_mask": mask_out}


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
