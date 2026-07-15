# The single point of contact with Python. Every other file goes through
# db_py() and py_try(); no file imports reticulate directly.
#
# Policy (see DESIGN.md section 7): running an analysis NEVER downloads or
# installs anything. diamondback uses a Python you already have, or it stops
# and tells you to run diamondback_install_python(). reticulate's uv-managed
# environment is genuinely convenient, but it will fetch a whole CPython
# interpreter plus NumPy and SciPy on first use, and a call to label_patches()
# is not consent to download 300 MB.

# Bump when the labelling kernel changes semantics. This, not the package
# version, invalidates caches -- a docstring fix should not throw away an hour
# of labelling. See DESIGN.md section 5.
ALGORITHM_VERSION <- "1"

# Minimums are what the kernels actually need: numpy for the bincount/unique
# paths, scipy for ndimage.label and distance_transform_edt.
PY_REQUIREMENTS <- c("numpy>=1.22", "scipy>=1.9")

.db_state <- new.env(parent = emptyenv())

#' Where consent to use a managed Python environment is recorded
#'
#' Written by [diamondback_install_python()] and by nothing else. Its presence
#' is what allows later sessions to use the managed environment without asking
#' again.
#' @noRd
db_consent_file <- function() {
  file.path(tools::R_user_dir("diamondback", "config"), "python-managed")
}

#' The configured Python preference
#'
#' `getOption("diamondback.python")`, then `DIAMONDBACK_PYTHON`, then "auto".
#' Recognised values: "auto", "managed", "system", or a path to a Python.
#' @noRd
db_python_option <- function() {
  opt <- getOption("diamondback.python", default = NULL)
  if (!is.null(opt) && nzchar(opt)) return(opt)
  ev <- Sys.getenv("DIAMONDBACK_PYTHON", "")
  if (nzchar(ev)) return(ev)
  "auto"
}

#' Find a Python that already has NumPy and SciPy
#'
#' Shells out to candidate interpreters rather than going through reticulate,
#' because initialising reticulate is itself what can trigger provisioning.
#' Returns a path, or NULL.
#' @noRd
db_probe_python <- function() {
  cands <- c(
    if (nzchar(Sys.getenv("VIRTUAL_ENV"))) {
      file.path(Sys.getenv("VIRTUAL_ENV"), if (.Platform$OS.type == "windows") "Scripts/python.exe" else "bin/python")
    },
    if (nzchar(Sys.getenv("CONDA_PREFIX"))) {
      file.path(Sys.getenv("CONDA_PREFIX"), if (.Platform$OS.type == "windows") "python.exe" else "bin/python")
    },
    unname(Sys.which("python3")),
    unname(Sys.which("python"))
  )
  cands <- unique(cands[nzchar(cands) & file.exists(cands)])

  for (p in cands) {
    ok <- tryCatch({
      out <- suppressWarnings(system2(p, c("-c", shQuote("import numpy, scipy")),
                                      stdout = TRUE, stderr = TRUE))
      is.null(attr(out, "status")) || identical(attr(out, "status"), 0L)
    }, error = function(e) FALSE)
    if (isTRUE(ok)) return(p)
  }
  NULL
}

#' Is a Python interpreter already running in this session?
#'
#' Its own function so the resolution logic can be tested: once any test has
#' initialised Python, every later call would otherwise short-circuit to
#' "existing" and the interesting branches would never run.
#' @noRd
db_py_initialized <- function() {
  isTRUE(tryCatch(reticulate::py_available(initialize = FALSE), error = function(e) FALSE))
}

#' Decide how Python will be obtained, without obtaining it
#'
#' Returns one of "existing", "user", "managed", "system", "none". Nothing here
#' initialises an interpreter or downloads anything.
#' @noRd
db_python_mode <- function() {
  if (db_py_initialized()) return("existing")

  mode <- db_python_option()

  # An explicit path is a user decision; honour it and never install over it.
  if (!mode %in% c("auto", "managed", "system")) {
    if (!file.exists(mode)) {
      cli::cli_abort(c(
        "The Python set in {.code diamondback.python} does not exist.",
        "x" = "{.path {mode}}"
      ), call = NULL)
    }
    Sys.setenv(RETICULATE_PYTHON = mode)
    return("user")
  }

  if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) return("user")

  # Consent, once given, persists.
  if (identical(mode, "managed") || file.exists(db_consent_file())) return("managed")

  p <- db_probe_python()
  if (!is.null(p)) {
    Sys.setenv(RETICULATE_PYTHON = p)
    return("system")
  }

  if (identical(mode, "system")) return("none")
  "none"
}

#' Load the diamondback Python backend
#'
#' Resolves and imports the Python backend. Most users never call this
#' directly; it runs on demand the first time a computation needs it. It is
#' exported so a session can be warmed up before a long run, and so that
#' [diamondback_check()] has something to report on.
#'
#' @section This never installs anything:
#' Running an analysis will not download or install software. diamondback looks,
#' in order, for:
#'
#' 1. a Python already initialised in this session;
#' 2. `RETICULATE_PYTHON`, or a path set in `options(diamondback.python = )`;
#' 3. the managed environment, if you have run [diamondback_install_python()];
#' 4. any Python on your `PATH` (or in an active virtualenv/conda environment)
#'    that already has NumPy and SciPy.
#'
#' If none of those works it stops with instructions. Downloading is only ever
#' done by [diamondback_install_python()], which asks first.
#'
#' @param quiet Suppress the startup message.
#' @return The Python module, invisibly, as a `python.builtin.module` object.
#' @seealso [diamondback_install_python()] to set up an environment,
#'   [diamondback_check()] for a full report.
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

  mode <- db_python_mode()
  if (identical(mode, "none")) db_python_missing_error()

  # py_require() is what enables reticulate's uv-managed environment, and with
  # it the download. It is called only when the user has opted in.
  if (identical(mode, "managed")) {
    try(reticulate::py_require(PY_REQUIREMENTS), silent = TRUE)
  }

  path <- system.file("python", package = "diamondback")
  if (!nzchar(path)) {
    cli::cli_abort("Could not locate the bundled Python module. Is diamondback installed correctly?")
  }

  mod <- tryCatch(
    reticulate::import_from_path("diamondback_core", path = path,
                                 delay_load = FALSE, convert = FALSE),
    error = function(e) db_python_setup_error(e, mode)
  )

  .db_state$py <- mod
  .db_state$mode <- mode
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

#' Is the Python backend usable?
#'
#' `TRUE` if diamondback can find a Python with NumPy and SciPy, `FALSE`
#' otherwise. Like everything else here except [diamondback_install_python()],
#' it never installs anything: it answers the question rather than fixing it.
#'
#' Useful for scripts, reports and examples that should degrade gracefully
#' rather than fail on a machine without the backend. The package's own examples
#' are guarded with it.
#'
#' @return A single logical.
#' @seealso [diamondback_check()] for a full report,
#'   [diamondback_install_python()] to set an environment up.
#' @export
#' @examples
#' if (diamondback_ready()) {
#'   label_patches(matrix(c(1, 0, 0, 1), 2, 2), quiet = TRUE)
#' }
diamondback_ready <- function() {
  isTRUE(tryCatch({
    diamondback_python(quiet = TRUE)
    TRUE
  }, error = function(e) FALSE))
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

# The no-data sentinel carried in int32 label blocks. It is R's NA_integer_ and
# the GeoTIFF NAflag, and it is one below the smallest usable patch ID.
DB_INT32_NA <- -2147483648

#' Convert a block of int32 label cells into an R integer vector
#'
#' reticulate does not convert numpy int32 identically on every platform: on
#' macOS/Linux it yields an R integer vector in which the INT_MIN sentinel has
#' already become NA_integer_, while on Windows it yields a double. The double
#' path still *works*, because as.integer(-2147483648) is NA -- but only by way
#' of an "NAs introduced by coercion to integer range" warning, once per block.
#' Correct by accident and noisy by design is not a contract worth keeping, so
#' the sentinel is mapped explicitly here and the platform difference stops at
#' this function.
#'
#' A value that is out of integer range and *not* the sentinel still warns,
#' because that would be a genuine bug rather than an expected encoding.
#' @noRd
db_int_rows <- function(rows) {
  v <- py_r(rows)
  if (is.integer(v)) return(as.vector(v))
  v <- as.vector(as.numeric(v))
  v[!is.na(v) & v == DB_INT32_NA] <- NA_real_
  as.integer(v)
}

db_python_missing_error <- function() {
  cli::cli_abort(
    c(
      "diamondback needs Python with NumPy and SciPy, and could not find one.",
      "x" = "No Python on your {.envvar PATH} has both installed.",
      "i" = "diamondback will not download anything on its own. Pick one:",
      "*" = "{.run diamondback_install_python()} -- set up a private environment \\
             for diamondback (downloads ~200 MB; asks first).",
      "*" = "Install them into the Python you already use: {.code pip install numpy scipy}.",
      "*" = "Point diamondback at an environment that has them: \\
             {.code options(diamondback.python = \"/path/to/python\")}."
    ),
    class = "diamondback_python_missing",
    call = NULL
  )
}

db_python_setup_error <- function(e, mode = "unknown") {
  msg <- conditionMessage(e)
  missing_mod <- grepl("No module named 'numpy'|No module named 'scipy'|ModuleNotFoundError", msg)

  where <- tryCatch(reticulate::py_config()$python, error = function(e) NA_character_)

  bullets <- c(
    "Could not load the diamondback Python backend.",
    "x" = if (missing_mod) {
      "NumPy and/or SciPy are missing from the Python being used."
    } else {
      "The Python environment could not be initialised."
    }
  )
  if (!is.na(where)) bullets <- c(bullets, "i" = "Python: {.path {where}}")
  if (missing_mod) {
    bullets <- c(
      bullets,
      "*" = "{.code pip install numpy scipy} into that environment, or",
      "*" = "{.run diamondback_install_python()} to set up a private one."
    )
  }
  bullets <- c(bullets, "i" = "{.run diamondback_check()} reports the full picture.",
               " " = "Python said: {msg}")

  cli::cli_abort(bullets, class = "diamondback_python_error", call = NULL)
}

#' Set up a private Python environment for diamondback
#'
#' The deliberate, opt-in installer. This is the **only** thing in the package
#' that downloads software; no analysis function will ever do it for you.
#'
#' It uses `reticulate::py_require()`, which resolves a `uv`-managed environment
#' containing a CPython interpreter, NumPy and SciPy. Your existing Python
#' installations, virtualenvs and conda environments are not touched or
#' modified. Consent is recorded in `tools::R_user_dir("diamondback", "config")`
#' so that later sessions use the environment without asking again.
#'
#' You do not need this if you already have a Python with NumPy and SciPy —
#' diamondback finds it on its own.
#'
#' @param prompt Ask for confirmation before downloading. Defaults to `TRUE` in
#'   interactive sessions. Set `FALSE` for scripted or CI setup, which is an
#'   explicit statement that the download is expected.
#' @param quiet Suppress progress reporting.
#' @return `TRUE` if the environment is ready, `FALSE` if you declined,
#'   invisibly.
#' @seealso [diamondback_check()], [diamondback_remove_python()]
#' @export
#' @examples
#' \dontrun{
#' diamondback_install_python()
#' }
diamondback_install_python <- function(prompt = interactive(), quiet = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    cli::cli_abort('The {.pkg reticulate} package is required. {.run install.packages("reticulate")}')
  }

  if (!quiet) {
    cli::cli_h2("Set up a Python environment for diamondback")
    cli::cli_text("This will download and cache, via {.pkg reticulate} and {.code uv}:")
    cli::cli_bullets(c(
      "*" = "a CPython interpreter, if one is not already cached",
      "*" = "{.val {PY_REQUIREMENTS}}"
    ))
    cli::cli_text("")
    cli::cli_bullets(c(
      "i" = "Roughly 200-350 MB on disk, under {.path {tools::R_user_dir('reticulate', 'cache')}}.",
      "i" = "Your existing Python, virtualenv and conda environments are {.strong not} modified.",
      "i" = "Remove it later with {.run diamondback_remove_python()}."
    ))
  }

  if (isTRUE(prompt)) {
    ans <- utils::menu(c("Yes, download and set it up", "No, cancel"),
                       title = "Download and set up this environment now?")
    if (!identical(ans, 1L)) {
      cli::cli_alert_info("Cancelled. Nothing was downloaded.")
      return(invisible(FALSE))
    }
  }

  # Record consent before resolving, so that an interrupted download does not
  # leave the user re-consenting on every attempt.
  f <- db_consent_file()
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "# diamondback: consent to use a reticulate/uv-managed Python environment.",
    "# Written by diamondback_install_python(). Delete this file, or run",
    "# diamondback_remove_python(), to revoke it.",
    paste("date:", as.character(Sys.time())),
    paste("requirements:", paste(PY_REQUIREMENTS, collapse = ", "))
  ), f)

  # Force a fresh resolution rather than reusing a handle from earlier in the
  # session, which may point at a different interpreter.
  .db_state$py <- NULL
  ok <- tryCatch({
    reticulate::py_require(PY_REQUIREMENTS)
    diamondback_python(quiet = TRUE)
    TRUE
  }, error = function(e) {
    unlink(f)
    cli::cli_abort(c(
      "Could not set up the Python environment.",
      "x" = conditionMessage(e),
      "i" = "Nothing was recorded; you can retry, or install NumPy and SciPy yourself."
    ), call = NULL)
  })

  if (!quiet && isTRUE(ok)) {
    v <- py_r(.db_state$py$versions())
    cli::cli_alert_success(
      "Ready: Python {v$python}, NumPy {v$numpy}, SciPy {v$scipy}."
    )
    cli::cli_alert_info("diamondback will use this environment from now on.")
  }
  invisible(TRUE)
}

#' Forget the managed Python environment
#'
#' Revokes the consent recorded by [diamondback_install_python()], so
#' diamondback stops using the managed environment. The cached files belong to
#' reticulate and are shared with other packages, so they are left alone; the
#' path to delete by hand is reported.
#'
#' @param quiet Suppress the report.
#' @return `TRUE` if consent was present and removed, invisibly.
#' @export
#' @examples
#' \dontrun{
#' diamondback_remove_python()
#' }
diamondback_remove_python <- function(quiet = FALSE) {
  f <- db_consent_file()
  had <- file.exists(f)
  if (had) unlink(f)
  .db_state$py <- NULL
  .db_state$mode <- NULL

  if (!quiet) {
    if (had) {
      cli::cli_alert_success("diamondback will no longer use the managed Python environment.")
      cli::cli_alert_info(
        "The cached files are reticulate's and may be shared with other packages. \\
         To reclaim the space, delete {.path {tools::R_user_dir('reticulate', 'cache')}}."
      )
    } else {
      cli::cli_alert_info("No managed environment was set up; nothing to remove.")
    }
  }
  invisible(had)
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
#' This never installs anything. If no suitable Python is found it says so and
#' tells you how to get one.
#'
#' @param verbose Print a formatted report. The value is returned invisibly
#'   either way, so this is also usable in scripts.
#' @return A list with elements `r`, `terra`, `reticulate`, `python`, `numpy`,
#'   `scipy`, `mode`, `memory_gb`, `test_passed` and `ok`, invisibly.
#' @seealso [diamondback_install_python()]
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
    mode = NA_character_,
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
    out$mode <- .db_state$mode %||% NA_character_
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
    if (!is.na(out$mode)) {
      cli::cli_text("{.emph Source:} {db_mode_label(out$mode)}")
    }
    cli::cli_h2("End-to-end test")
    if (out$ok) {
      cli::cli_alert_success("Labelling a small test raster succeeded. diamondback is ready.")
    } else {
      cli::cli_alert_danger("Test labelling failed.")
      if (!is.null(out$error)) cli::cli_text("{.emph Error:} {out$error}")
    }
  }

  invisible(out)
}

db_mode_label <- function(mode) {
  switch(mode,
    existing = "a Python already active in this session",
    user = "your RETICULATE_PYTHON / diamondback.python setting",
    managed = "diamondback's managed environment (diamondback_install_python)",
    system = "a Python found on your PATH that already had NumPy and SciPy",
    mode
  )
}
