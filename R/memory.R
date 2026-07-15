# Memory estimation. SciPy's label() needs the array in RAM, so the only
# honest thing to do is to say what that will cost before allocating it.

# Peak bytes per cell for each stage, from DESIGN.md section 10.
#   label   : code uint8 (1) + boolean foreground (1) + int32 labels (4)
#   metrics : code (1) + labels (4), plus chunk-sized temporaries
#   core    : code (1) + labels (4) + float64 EDT distances (8), minus overlap
BYTES_PER_CELL <- c(label = 6, metrics = 5, core = 13)

#' Detect available system memory in bytes
#'
#' Falls back to `NA` rather than guessing, because a wrong guess here turns a
#' helpful check into a spurious blocker.
#' @noRd
db_available_memory <- function() {
  sys <- Sys.info()[["sysname"]]
  val <- tryCatch({
    if (identical(sys, "Darwin")) {
      # Free + inactive + speculative pages: memory the OS can hand out now.
      page <- as.numeric(system2("sysctl", c("-n", "hw.pagesize"), stdout = TRUE))
      vm <- system2("vm_stat", stdout = TRUE)
      get <- function(key) {
        line <- grep(key, vm, value = TRUE, fixed = TRUE)
        if (!length(line)) return(0)
        as.numeric(gsub("[^0-9]", "", line[1]))
      }
      (get("Pages free:") + get("Pages inactive:") + get("Pages speculative:")) * page
    } else if (identical(sys, "Linux")) {
      mi <- readLines("/proc/meminfo", warn = FALSE)
      line <- grep("^MemAvailable:", mi, value = TRUE)
      if (!length(line)) line <- grep("^MemFree:", mi, value = TRUE)
      as.numeric(gsub("[^0-9]", "", line[1])) * 1024
    } else if (identical(sys, "Windows")) {
      # utils::memoryLimit() exists only on Windows, so it is fetched
      # dynamically: naming it directly makes R CMD check flag it as missing
      # on every other platform.
      ml <- tryCatch({
        f <- get("memoryLimit", envir = asNamespace("utils"), inherits = FALSE)
        f()
      }, error = function(e) NA_real_)
      if (is.numeric(ml) && length(ml) == 1 && !is.na(ml) && ml > 0 && is.finite(ml)) {
        ml * 1024^2
      } else {
        out <- system2("wmic", c("OS", "get", "FreePhysicalMemory", "/Value"), stdout = TRUE)
        line <- grep("FreePhysicalMemory", out, value = TRUE)
        as.numeric(gsub("[^0-9]", "", line[1])) * 1024
      }
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_, warning = function(w) NA_real_)

  if (!is.numeric(val) || length(val) != 1 || is.na(val) || val <= 0) NA_real_ else val
}

db_fmt_bytes <- function(x) {
  if (is.na(x)) return("unknown")
  if (x >= 1024^3) return(sprintf("%.1f GB", x / 1024^3))
  if (x >= 1024^2) return(sprintf("%.0f MB", x / 1024^2))
  sprintf("%.0f KB", x / 1024)
}

#' Estimate peak memory for a diamondback stage
#'
#' @param ncell Number of cells that will be held in memory.
#' @param stage One of "label", "metrics", "core".
#' @return Estimated peak bytes.
#' @noRd
db_estimate_memory <- function(ncell, stage = "label") {
  stage <- match.arg(stage, names(BYTES_PER_CELL))
  as.numeric(ncell) * BYTES_PER_CELL[[stage]]
}

#' Stop before an allocation that will not fit
#'
#' Errors *before* anything is allocated, naming the estimate, the ceiling, and
#' the ways out. Never guesses: if available memory cannot be detected, this
#' warns at an absolute threshold rather than blocking work.
#'
#' @param ncell Cells to be held in memory.
#' @param stage Stage name, used for the per-cell cost and the message.
#' @param max_memory_frac Fraction of available RAM that may be used.
#' @param memory_limit Explicit ceiling in bytes, overriding detection.
#' @noRd
db_check_memory <- function(ncell, stage = "label", max_memory_frac = 0.6,
                            memory_limit = NULL, quiet = FALSE) {
  need <- db_estimate_memory(ncell, stage)

  avail <- if (!is.null(memory_limit)) as.numeric(memory_limit) else db_available_memory()
  ceiling_bytes <- if (is.na(avail)) NA_real_ else avail * max_memory_frac

  if (!quiet) {
    cli::cli_alert_info(
      "Estimated peak memory for {stage}: {db_fmt_bytes(need)} \\
       ({.val {format(ncell, big.mark = ',')}} cells)."
    )
  }

  if (is.na(ceiling_bytes)) {
    if (need > 8 * 1024^3) {
      cli::cli_warn(c(
        "Could not detect available memory, and this operation needs about {db_fmt_bytes(need)}.",
        "i" = "If it fails, crop or mask to a smaller domain, or set {.arg memory_limit}."
      ))
    }
    return(invisible(need))
  }

  if (need > ceiling_bytes) {
    cli::cli_abort(c(
      "This operation needs about {db_fmt_bytes(need)}, which exceeds the safe limit of {db_fmt_bytes(ceiling_bytes)}.",
      "x" = "{stage} needs roughly {BYTES_PER_CELL[[stage]]} bytes per cell for \\
             {.val {format(ncell, big.mark = ',')}} cells.",
      "i" = "Detected {db_fmt_bytes(avail)} available; the limit is {max_memory_frac} of that.",
      "*" = "Crop or mask to a smaller analysis domain (the usual answer).",
      "*" = "Raise {.arg max_memory_frac} if you know the estimate is pessimistic.",
      "*" = "Set {.arg memory_limit} to override memory detection entirely."
    ), class = "diamondback_memory_error", call = NULL)
  }

  if (need > ceiling_bytes * 0.5 && !quiet) {
    cli::cli_alert_warning(
      "This will use a large share of available memory ({db_fmt_bytes(need)} of {db_fmt_bytes(avail)})."
    )
  }
  invisible(need)
}

#' Report estimated memory cost of analysing a raster
#'
#' Prints what each diamondback stage is expected to cost for a given raster,
#' without doing any work. Useful before committing to a long run.
#'
#' @param x A `SpatRaster`, raster filename, or matrix.
#' @return A named numeric vector of estimated peak bytes per stage, invisibly.
#' @export
#' @examples
#' r <- terra::rast(nrows = 100, ncols = 100)
#' terra::values(r) <- rep(c(0, 1), terra::ncell(r) / 2)
#' db_memory_report(r)
db_memory_report <- function(x) {
  r <- db_as_rast(x)
  n <- terra::ncell(r)
  est <- vapply(names(BYTES_PER_CELL), function(s) db_estimate_memory(n, s), numeric(1))
  avail <- db_available_memory()

  cli::cli_h2("Memory estimate")
  cli::cli_text("Raster: {.val {terra::nrow(r)}} x {.val {terra::ncol(r)}} = \\
                 {.val {format(n, big.mark = ',')}} cells")
  cli::cli_dl(c(
    "label_patches()" = db_fmt_bytes(est[["label"]]),
    "patch_metrics()" = db_fmt_bytes(est[["metrics"]]),
    "patch_core_area()" = db_fmt_bytes(est[["core"]]),
    "Available RAM" = db_fmt_bytes(avail)
  ))
  invisible(est)
}
