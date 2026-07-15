# Memory estimation. SciPy's label() needs the array in RAM, so the only
# honest thing to do is to say what that will cost before allocating it.

# Peak bytes per cell for each stage, from DESIGN.md section 10. These assume a
# uint8 code array; runs with more than 253 classes use uint16, and
# `code_bytes` accounts for the extra byte.
#   label   : code (1) + boolean foreground (1) + int32 labels (4)
#   metrics : code (1) + labels (4), plus chunk-sized temporaries
#   core    : code (1) + labels (4), plus a bounded row strip. The full-array
#             float64 distance grid that used to make this 13 bytes/cell is gone;
#             the transform is tiled, so the distances never exist all at once.
BYTES_PER_CELL <- c(label = 6, metrics = 5, core = 6)

#' Total physical memory in bytes
#'
#' The denominator for "could this ever fit". Distinct from
#' `db_available_memory()`, which answers "is it free *right now*" -- a much
#' smaller and much more volatile number, and the wrong one to refuse work over.
#' @noRd
db_physical_memory <- function() {
  sys <- Sys.info()[["sysname"]]
  val <- tryCatch({
    if (identical(sys, "Darwin")) {
      as.numeric(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE))
    } else if (identical(sys, "Linux")) {
      mi <- readLines("/proc/meminfo", warn = FALSE)
      line <- grep("^MemTotal:", mi, value = TRUE)
      if (!length(line)) return(NA_real_)
      as.numeric(gsub("[^0-9]", "", line[1])) * 1024
    } else if (identical(sys, "Windows")) {
      out <- system2("wmic", c("ComputerSystem", "get", "TotalPhysicalMemory", "/Value"),
                     stdout = TRUE)
      line <- grep("TotalPhysicalMemory", out, value = TRUE)
      if (!length(line)) return(NA_real_)
      as.numeric(gsub("[^0-9]", "", line[1]))
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_, warning = function(w) NA_real_)

  if (!is.numeric(val) || length(val) != 1 || is.na(val) || val <= 0) NA_real_ else val
}

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
db_estimate_memory <- function(ncell, stage = "label", code_bytes = 1) {
  stage <- match.arg(stage, names(BYTES_PER_CELL))
  as.numeric(ncell) * (BYTES_PER_CELL[[stage]] + (code_bytes - 1))
}

#' Refuse only what cannot work; warn about the rest
#'
#' Memory availability cannot be known, so this does not pretend to. It answers
#' two different questions with two different numbers:
#'
#' * **Can this fit at all?** Compared against *physical* RAM. Only this may
#'   refuse, because only this is structural: no amount of closing browser tabs
#'   makes a 40 GB array fit in 8 GB.
#' * **Will it be comfortable?** Compared against *currently free* memory. This
#'   may only warn. The OS reclaims, compresses and pages; free memory is a
#'   snapshot, not a budget.
#'
#' The first design refused whenever the estimate exceeded a fraction of free
#' memory, and had it backwards: on macOS it read 1.7 GB free on an 8 GB machine
#' and would have blocked a 4 GB job that runs there routinely. A false refusal
#' costs the whole analysis; a warning costs only attention.
#'
#' @param ncell Cells to be held in memory.
#' @param stage Stage name, used for the per-cell cost and the message.
#' @param max_memory_frac Fraction of *physical* RAM above which to hard-stop.
#' @param memory_limit Explicit ceiling in bytes, replacing physical RAM.
#' @noRd
db_check_memory <- function(ncell, stage = "label", max_memory_frac = 0.9,
                            memory_limit = NULL, quiet = FALSE, code_bytes = 1) {
  need <- db_estimate_memory(ncell, stage, code_bytes)

  physical <- if (!is.null(memory_limit)) as.numeric(memory_limit) else db_physical_memory()
  free <- db_available_memory()
  hard <- if (is.na(physical)) NA_real_ else physical * max_memory_frac

  if (!quiet) {
    cli::cli_alert_info(
      "Estimated peak memory for {stage}: {db_fmt_bytes(need)} \\
       ({.val {format(ncell, big.mark = ',')}} cells)."
    )
  }

  # The only refusal: it cannot fit even with the machine to itself.
  if (!is.na(hard) && need > hard) {
    cli::cli_abort(c(
      "This operation needs about {db_fmt_bytes(need)}, more than this machine has.",
      "x" = "{stage} needs roughly {BYTES_PER_CELL[[stage]] + code_bytes - 1} bytes per \\
             cell for {.val {format(ncell, big.mark = ',')}} cells.",
      "i" = "Physical memory is {db_fmt_bytes(physical)}; the ceiling is \\
             {max_memory_frac} of that ({db_fmt_bytes(hard)}).",
      "*" = "Crop or mask to a smaller analysis domain (the usual answer).",
      "*" = "Raise {.arg max_memory_frac}, or set {.arg memory_limit}, if you know better.",
      "i" = "This is a structural limit, not a guess about what is free right now."
    ), class = "diamondback_memory_error", call = NULL)
  }

  if (is.na(physical)) {
    if (need > 8 * 1024^3) {
      cli::cli_warn(c(
        "Could not detect physical memory, and this needs about {db_fmt_bytes(need)}.",
        "i" = "If it fails, crop or mask to a smaller domain, or set {.arg memory_limit}."
      ))
    }
    return(invisible(need))
  }

  if (!quiet) {
    if (!is.na(free) && need > free) {
      # Not a problem, just slow: the OS will page. Worth saying, because a run
      # that suddenly takes ten times longer is otherwise a mystery.
      cli::cli_alert_warning(
        "This needs {db_fmt_bytes(need)} but only {db_fmt_bytes(free)} is free right \\
         now, so expect paging and a slower run. Closing other applications would help."
      )
    } else if (need > physical * 0.5) {
      cli::cli_alert_warning(
        "This will use a large share of this machine's memory \\
         ({db_fmt_bytes(need)} of {db_fmt_bytes(physical)})."
      )
    }
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
