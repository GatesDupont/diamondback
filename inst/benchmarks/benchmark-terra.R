# Benchmark: diamondback vs terra::patches()
#
# Deliberately NOT part of the test suite. Timings depend on hardware, disk
# and what else the machine is doing, so asserting on them would produce a
# suite that fails for reasons that have nothing to do with the code.
#
# Run with:
#   Rscript inst/benchmarks/benchmark-terra.R           # default sizes
#   Rscript inst/benchmarks/benchmark-terra.R 4000      # up to 4000x4000
#
# Correctness is checked as it goes: any size where the two disagree on the
# patch count is reported loudly, because a fast wrong answer is worthless.

suppressMessages({
  library(terra)
  library(diamondback)
})

args <- commandArgs(trailingOnly = TRUE)
max_side <- if (length(args)) as.numeric(args[[1]]) else 2000
sides <- c(250, 500, 1000, 2000, 4000, 8000)
sides <- sides[sides <= max_side]
density <- 0.4

make_raster <- function(side, seed = 1) {
  r <- rast(nrows = side, ncols = side, xmin = 0, xmax = side * 30,
            ymin = 0, ymax = side * 30, crs = "EPSG:5070")
  set.seed(seed)
  # Smoothed noise gives patch structure like real land cover, rather than the
  # salt-and-pepper that pure rbinom produces.
  values(r) <- rbinom(ncell(r), 1, density)
  r <- focal(r, w = 3, fun = "mean", na.rm = TRUE)
  r <- ifel(r > 0.45, 1, 0)
  r
}

timeit <- function(expr) {
  gc(verbose = FALSE)
  t <- system.time(force(expr))
  as.numeric(t[["elapsed"]])
}

cat("\ndiamondback vs terra::patches()\n")
cat(sprintf("%-12s %12s %12s %12s %10s %10s %8s\n",
            "size", "cells", "diamondback", "terra", "speedup", "patches", "agree"))
cat(strrep("-", 82), "\n")

results <- list()

for (side in sides) {
  r <- make_raster(side)
  n <- ncell(r)

  # diamondback: labelling only, to compare like with like.
  t_db <- timeit(db <- label_patches(r, directions = 8, quiet = TRUE))
  n_db <- db$metadata$n_patches

  # terra wants NA background rather than 0.
  rt <- r
  rt[rt == 0] <- NA
  t_tr <- timeit(tr <- patches(rt, directions = 8, zeroAsNA = FALSE))
  n_tr <- length(unique(na.omit(values(tr, mat = FALSE))))

  agree <- identical(as.integer(n_db), as.integer(n_tr))
  cat(sprintf("%-12s %12s %12.2f %12.2f %9.1fx %10s %8s\n",
              sprintf("%dx%d", side, side),
              format(n, big.mark = ","),
              t_db, t_tr, t_tr / t_db,
              format(n_db, big.mark = ","),
              if (agree) "yes" else "NO!"))

  if (!agree) {
    cat(sprintf("   !! MISMATCH: diamondback %d vs terra %d\n", n_db, n_tr))
  }

  results[[length(results) + 1L]] <- data.frame(
    side = side, cells = n, diamondback = t_db, terra = t_tr,
    speedup = t_tr / t_db, patches = n_db, agree = agree
  )
  rm(r, rt, tr, db); gc(verbose = FALSE)
}

cat("\n")

# --- where the time actually goes ------------------------------------------
side <- max(sides)
r <- make_raster(side)
cat(sprintf("Stage breakdown at %dx%d (%s cells)\n", side, side,
            format(ncell(r), big.mark = ",")))
t_label <- timeit(lab <- label_patches(r, directions = 8, quiet = TRUE))
t_metrics <- timeit(met <- patch_metrics(lab, quiet = TRUE))
t_core <- timeit(patch_core_area(met, depth = 60, quiet = TRUE))
t_cmp <- timeit(compare_patches(lab, lab, quiet = TRUE))
cat(sprintf("  label_patches()   %6.2f s\n", t_label))
cat(sprintf("  patch_metrics()   %6.2f s\n", t_metrics))
cat(sprintf("  patch_core_area() %6.2f s\n", t_core))
cat(sprintf("  compare_patches() %6.2f s\n", t_cmp))

cat("\nNotes\n")
cat("  * terra::patches() timings include its own IO; both run on an in-memory raster.\n")
cat("  * diamondback's first call in a session also pays Python startup (~1-2 s),\n")
cat("    which is excluded here by warming up above but matters for one-shot scripts.\n")

invisible(do.call(rbind, results))
