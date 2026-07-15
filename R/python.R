# The single point of contact with Python. Every other file goes through
# db_py() and py_try(); no file imports reticulate directly.

# Bump when the labelling kernel changes semantics. This, not the package
# version, invalidates caches -- a docstring fix should not throw away an hour
# of labelling. See DESIGN.md section 5.
ALGORITHM_VERSION <- "1"

# Minimums are what the kernels actually need: numpy for the bincount/unique
# paths, scipy for ndimage.label and distance_transform_edt.
PY_REQUIREMENTS <- c("numpy>=1.22", "scipy>=1.9")

.db_state <- new.env(parent = emptyenv())

#' Load the diamondback Python backend
#'
#' Imports the bundled Python module, declaring the NumPy and SciPy
#' requirements on first use. Nothing is installed into your existing Python
#' environments and nothing happens when the package is loaded: the backend is
#' resolved the first time a computation actually needs it.
#'
#' Most users never call this directly. It is exported so that a session can be
#' warmed up ahead of a long run, and so [diamondback_check()] has something to
#' report on.
#'
#' @param quiet Suppress the startup message.
#' @return The Python module, invisibly, as a `python.builtin.module` object.
#' @seealso [diamondback_check()] for a full environment report.
#' @export
#' @examples
#' \dontrun{
#' diamondback_python()
#' }
diamondback_python <- function(quiet = TRUE) {
  if (!is.null(.db_state$py)) {
    return(invisible(.db_state$py))
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg reticulate} package is required but not installed.",
      "i" = 'Install it with {.run install.packages("reticulate")}.'
    ))
  }

  # py_require() (reticulate >= 1.41) declares dependencies and lets reticulate
  # resolve an ephemeral uv-managed environment. If the user has pointed
  # RETICULATE_PYTHON at their own environment, it steps aside and we simply
  # check what is there.
  if (!nzchar(Sys.getenv("RETICULATE_PYTHON", unset = ""))) {
    try(reticulate::py_require(PY_REQUIREMENTS), silent = TRUE)
  }

  path <- system.file("python", package = "diamondback")
  if (!nzchar(path)) {
    cli::cli_abort("Could not locate the bundled Python module. Is diamondback installed correctly?")
  }

  # convert = FALSE is load-bearing, not a preference. With reticulate's default
  # every return value is converted to an R object, which would pull each label
  # array back into R and defeat the entire point of doing the work in NumPy.
  # Big arrays therefore stay as Python handles; small per-patch results are
  # converted explicitly with py_r().
  mod <- tryCatch(
    reticulate::import_from_path("diamondback_core", path = path,
                                 delay_load = FALSE, convert = FALSE),
    error = function(e) db_python_setup_error(e)
  )

  .db_state$py <- mod
  if (!quiet) {
    v <- py_r(py_try(mod$versions(), "reading backend versions"))
    cli::cli_alert_success(
      "diamondback backend ready (Python {v$python}, NumPy {v$numpy}, SciPy {v$scipy})."
    )
  }
  invisible(mod)
}

# Internal shorthand.
db_py <- function() {
  if (is.null(.db_state$py)) diamondback_python(quiet = TRUE)
  .db_state$py
}

# The backend module is imported with convert = FALSE, so results come back as
# Python handles. These are the two ways to cross back into R.

# py_r(): convert a small Python result (a dict of per-patch vectors, a scalar)
# into R. Never call this on a full-raster array.
py_r <- function(x) {
  if (!inherits(x, "python.builtin.object")) return(x)
  reticulate::py_to_r(x)
}

# py_get(): pull one item out of a Python dict, keeping it in Python. This is
# how a big array is extracted without converting it.
py_get <- function(x, key) reticulate::py_get_item(x, key)

# py_num(): a numeric vector from one entry of a Python dict.
py_num <- function(x, key) as.numeric(py_r(py_get(x, key)))

db_python_setup_error <- function(e) {
  msg <- conditionMessage(e)
  hint <- if (grepl("No module named 'numpy'|No module named 'scipy'|ModuleNotFoundError", msg)) {
    "NumPy and/or SciPy are missing from the Python environment being used."
  } else {
    "The Python environment could not be initialised."
  }
  cli::cli_abort(
    c(
      "Could not load the diamondback Python backend.",
      "x" = hint,
      "i" = "Run {.run diamondback_check()} for a full diagnostic.",
      "i" = "diamondback needs Python with {.val {PY_REQUIREMENTS}}.",
      " " = "Python said: {msg}"
    ),
    class = "diamondback_python_error",
    call = NULL
  )
}

#' Translate Python exceptions into useful R errors
#'
#' Wraps a reticulate call so that a Python traceback becomes an R condition
#' naming the diamondback operation that failed. `MemoryError` is special-cased,
#' because a bare traceback is unhelpful when the real answer is "the raster did
#' not fit".
#'
#' @param expr Expression calling into the Python backend.
#' @param what Human-readable description of the operation, used in the message.
#' @return The value of `expr`.
#' @noRd
py_try <- function(expr, what = "a Python computation") {
  withCallingHandlers(
    tryCatch(
      expr,
      error = function(e) {
        msg <- conditionMessage(e)
        etype <- sub("^([A-Za-z_.]*Error|[A-Za-z_.]*Exception).*", "\\1",
                     strsplit(msg, "\n")[[1]][1])

        if (grepl("MemoryError|Unable to allocate|bad_alloc", msg)) {
          cli::cli_abort(c(
            "Ran out of memory during {what}.",
            "x" = "Python could not allocate the array it needed.",
            "i" = "Crop or mask the raster to a smaller analysis domain, or run on a machine with more RAM.",
            "i" = "{.fn db_memory_report} shows what a raster is estimated to cost.",
            " " = "Python said: {strsplit(msg, '\n')[[1]][1]}"
          ), class = c("diamondback_memory_error", "diamondback_python_error"), call = NULL)
        }

        cli::cli_abort(c(
          "Python backend failed during {what}.",
          "x" = "{etype}",
          "i" = "This is a bug in diamondback unless the inputs were unusual; please report it.",
          " " = "Python said: {msg}"
        ), class = "diamondback_python_error", call = NULL)
      }
    ),
    warning = function(w) {
      # NumPy emits RuntimeWarnings that are noise in this context.
      if (grepl("invalid value encountered|divide by zero", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

#' Report on the diamondback computing environment
#'
#' Checks everything the package needs and runs a small labelling operation
#' end-to-end, so that a green result means the whole R-to-Python path works,
#' not merely that the pieces are installed.
#'
#' Call this first when something is not working, and when setting up on a new
#' machine.
#'
#' @param verbose Print a formatted report. The value is returned invisibly
#'   either way, so this is also usable in scripts.
#' @return A list with elements `r`, `terra`, `reticulate`, `python`, `numpy`,
#'   `scipy`, `memory_gb`, `test_passed` and `ok`, invisibly.
#' @export
#' @examples
#' \dontrun{
#' diamondback_check()
#' }
diamondback_check <- function(verbose = TRUE) {
  out <- list(
    r = paste(R.version$major, R.version$minor, sep = "."),
    terra = tryCatch(as.character(utils::packageVersion("terra")), error = function(e) NA_character_),
    reticulate = tryCatch(as.character(utils::packageVersion("reticulate")), error = function(e) NA_character_),
    python = NA_character_,
    python_path = NA_character_,
    numpy = NA_character_,
    scipy = NA_character_,
    memory_gb = tryCatch(db_available_memory() / 1024^3, error = function(e) NA_real_),
    test_passed = FALSE,
    error = NULL
  )

  err <- NULL
  ok <- tryCatch({
    mod <- diamondback_python(quiet = TRUE)
    v <- py_r(mod$versions())
    out$python <- v$python
    out$numpy <- v$numpy
    out$scipy <- v$scipy
    out$python_path <- tryCatch(reticulate::py_config()$python, error = function(e) NA_character_)

    # An end-to-end labelling of a known 3x3 case: two 4-connected patches.
    m <- matrix(c(1, 0, 1,
                  0, 0, 0,
                  1, 1, 0), nrow = 3, byrow = TRUE)
    res <- label_patches(m, directions = 4, quiet = TRUE)
    out$test_passed <- isTRUE(nrow(res$metrics) == 3L)
    TRUE
  }, error = function(e) {
    err <<- conditionMessage(e)
    FALSE
  })

  out$error <- err
  out$ok <- isTRUE(ok) && isTRUE(out$test_passed)

  if (verbose) {
    cli::cli_h1("diamondback environment")
    cli::cli_dl(c(
      "R" = out$r,
      "terra" = out$terra,
      "reticulate" = out$reticulate,
      "Python" = if (is.na(out$python)) "{.strong not available}" else out$python,
      "NumPy" = if (is.na(out$numpy)) "{.strong not available}" else out$numpy,
      "SciPy" = if (is.na(out$scipy)) "{.strong not available}" else out$scipy,
      "Available RAM" = if (is.na(out$memory_gb)) "unknown" else sprintf("%.1f GB", out$memory_gb)
    ))
    if (!is.na(out$python_path)) {
      cli::cli_text("{.emph Python at} {.path {out$python_path}}")
    }
    cli::cli_h2("End-to-end test")
    if (out$ok) {
      cli::cli_alert_success("Labelling a small test raster succeeded. diamondback is ready.")
    } else {
      cli::cli_alert_danger("Test labelling failed.")
      if (!is.null(out$error)) cli::cli_text("{.emph Error:} {out$error}")
      cli::cli_alert_info("diamondback needs Python with {.val {PY_REQUIREMENTS}}.")
      cli::cli_alert_info(
        "By default reticulate resolves these automatically. To use your own \\
         environment instead, set {.envvar RETICULATE_PYTHON} before loading the package."
      )
    }
  }

  invisible(out)
}
