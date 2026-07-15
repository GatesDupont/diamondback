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

# Must equal INTERFACE_VERSION in inst/python/diamondback_core.py. Bump both
# together whenever the R-facing surface of that module changes: function names,
# signatures, or the keys of a returned dict.
PY_INTERFACE_VERSION <- "2"

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

#' Interpreter paths inside an environment prefix
#' @noRd
db_python_bin <- function(prefix) {
  if (.Platform$OS.type == "windows") {
    c(file.path(prefix, "python.exe"), file.path(prefix, "Scripts", "python.exe"))
  } else {
    file.path(prefix, "bin", "python")
  }
}

#' Plausible conda/mamba installation roots
#'
#' `reticulate::conda_list()` is the obvious source and cannot be the only one:
#' on a Homebrew miniforge layout it returns a single bogus entry named "envs"
#' pointing at a python that does not exist, and never sees the real
#' environments. So the roots are also searched directly. Everything found is
#' verified by actually importing NumPy and SciPy, so a wrong guess here costs
#' nothing but a moment.
#' @noRd
db_conda_roots <- function() {
  roots <- character()

  ce <- Sys.getenv("CONDA_EXE", "")            # <root>/bin/conda
  if (nzchar(ce)) roots <- c(roots, dirname(dirname(ce)))

  cp <- Sys.getenv("CONDA_PREFIX", "")
  if (nzchar(cp)) {
    roots <- c(roots, cp)                       # base, or an active env
    if (identical(basename(dirname(cp)), "envs")) {
      roots <- c(roots, dirname(dirname(cp)))   # the root above an active env
    }
  }

  roots <- c(
    roots,
    path.expand(c("~/miniforge3", "~/mambaforge", "~/miniconda3", "~/anaconda3",
                  "~/opt/anaconda3", "~/opt/miniconda3")),
    "/opt/homebrew/Caskroom/miniforge/base",
    "/opt/conda"
  )
  unique(roots[nzchar(roots) & dir.exists(roots)])
}

#' Every Python interpreter we can find, each with the name you would call it
#'
#' Ordered by how strong a signal the user gave: an activated environment first,
#' then whatever `python` means on the PATH, then conda environments, then
#' virtualenvs. Environment lists are sorted so that auto-detection is
#' deterministic rather than filesystem-order dependent.
#'
#' Nothing here initialises reticulate: doing so is what can trigger
#' provisioning, which must never happen by accident.
#' @noRd
db_python_candidates <- function() {
  cand <- list()
  push <- function(name, prefix = NULL, python = NULL) {
    ps <- if (is.null(python)) db_python_bin(prefix) else python
    ps <- ps[nzchar(ps) & file.exists(ps)]
    if (length(ps)) {
      cand[[length(cand) + 1L]] <<- list(name = name, python = unname(ps[1]))
    }
  }

  # An activated environment is an explicit statement of intent.
  ve <- Sys.getenv("VIRTUAL_ENV", "")
  if (nzchar(ve)) push(basename(ve), prefix = ve)
  cp <- Sys.getenv("CONDA_PREFIX", "")
  if (nzchar(cp)) push(basename(cp), prefix = cp)

  # Whatever "python" already means here.
  for (exe in c("python3", "python")) {
    w <- unname(Sys.which(exe))
    if (nzchar(w)) push(exe, python = w)
  }

  # conda: the base install, then each environment under it.
  for (root in db_conda_roots()) {
    push(basename(root), prefix = root)
    envs <- tryCatch(
      sort(list.dirs(file.path(root, "envs"), recursive = FALSE, full.names = TRUE)),
      error = function(e) character()
    )
    for (e in envs) push(basename(e), prefix = e)
  }

  # Anything reticulate knows about that actually exists.
  cl <- tryCatch(reticulate::conda_list(), error = function(e) NULL)
  if (is.data.frame(cl) && nrow(cl)) {
    for (i in seq_len(nrow(cl))) push(cl$name[i], python = cl$python[i])
  }

  # virtualenvs
  wh <- Sys.getenv("WORKON_HOME", "")
  if (!nzchar(wh)) wh <- path.expand("~/.virtualenvs")
  if (dir.exists(wh)) {
    for (e in sort(list.dirs(wh, recursive = FALSE, full.names = TRUE))) {
      push(basename(e), prefix = e)
    }
  }

  if (!length(cand)) return(cand)
  paths <- vapply(cand, `[[`, character(1), "python")
  cand[!duplicated(paths)]
}

#' Does this interpreter have what the kernels need?
#'
#' Shells out rather than going through reticulate, which can only be pointed at
#' one interpreter per session and whose initialisation is what may provision.
#' @noRd
db_has_deps <- function(python) {
  isTRUE(tryCatch({
    out <- suppressWarnings(system2(python, c("-c", shQuote("import numpy, scipy")),
                                    stdout = TRUE, stderr = TRUE))
    is.null(attr(out, "status")) || identical(attr(out, "status"), 0L)
  }, error = function(e) FALSE))
}

#' Find a Python that already has NumPy and SciPy
#'
#' Returns a path, or NULL. Any environment with NumPy and SciPy serves equally
#' well -- diamondback only uses them for its own kernels, not for your analysis
#' -- so taking the first that qualifies is not a compromise, and
#' [diamondback_check()] always reports which one was taken.
#' @noRd
db_probe_python <- function() {
  for (cd in db_python_candidates()) {
    if (db_has_deps(cd$python)) return(cd$python)
  }
  NULL
}

#' Resolve an environment *name* to an interpreter path
#'
#' So that `options(diamondback.python = "geo")` works and nobody has to go
#' hunting for `.../envs/geo/bin/python`. Returns a path, or NULL if no
#' environment goes by that name.
#' @noRd
db_resolve_env_name <- function(name) {
  for (cd in db_python_candidates()) {
    if (identical(cd$name, name)) return(cd$python)
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

  # Anything that is not a keyword is a user decision: a path, or the name of an
  # environment. Honour it, and never install over it.
  if (!mode %in% c("auto", "managed", "system")) {
    py <- if (file.exists(mode)) mode else db_resolve_env_name(mode)

    if (is.null(py)) {
      looks_like_path <- grepl("[/\\\\]", mode) || grepl("\\.exe$", mode)
      if (looks_like_path) {
        cli::cli_abort(c(
          "The Python set in {.code diamondback.python} does not exist.",
          "x" = "{.path {mode}}"
        ), call = NULL)
      }
      db_env_name_error(mode)
    }

    if (!db_has_deps(py)) {
      cli::cli_abort(c(
        "The Python set in {.code diamondback.python} is missing NumPy or SciPy.",
        "x" = "{.path {py}}",
        "i" = "Install them there, e.g. {.code {py} -m pip install numpy scipy}, \\
               or point at a different environment."
      ), class = "diamondback_python_missing", call = NULL)
    }

    Sys.setenv(RETICULATE_PYTHON = py)
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
#' 2. `RETICULATE_PYTHON`, or `options(diamondback.python = )`;
#' 3. the managed environment, if you have run [diamondback_install_python()];
#' 4. any environment it can find that already has NumPy and SciPy --- an active
#'    virtualenv or conda environment, `python` on your `PATH`, any conda
#'    environment, any virtualenv under `WORKON_HOME`.
#'
#' If none of those works it stops, **lists the environments it found**, and
#' tells you how to proceed. Downloading is only ever done by
#' [diamondback_install_python()], which asks first.
#'
#' @section Naming an environment:
#' `options(diamondback.python = )` takes an environment **name**, not just a
#' path, so you rarely need to know where anything lives:
#'
#' ```r
#' options(diamondback.python = "geo")     # a conda env or virtualenv called "geo"
#' options(diamondback.python = "python3") # whatever python3 means on your PATH
#' options(diamondback.python = "/opt/envs/geo/bin/python")  # still fine
#' ```
#'
#' A name is resolved against active environments, your `PATH`, conda
#' environments (found from `CONDA_EXE`, `CONDA_PREFIX`, and the usual install
#' roots) and virtualenvs. Whatever it resolves to is verified to have NumPy and
#' SciPy before it is used, and if the name matches nothing you get the list of
#' names that do exist. `"auto"`, `"managed"` and `"system"` are reserved.
#'
#' If you have exactly one environment with NumPy and SciPy, you can usually
#' skip this entirely: it will be found.
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

  # Always, not just when reporting: this is the handshake, not a nicety.
  v <- py_r(py_try(mod$versions(), "reading backend versions"))
  db_check_interface(v, path)

  .db_state$py <- mod
  .db_state$mode <- mode
  if (!quiet) {
    cli::cli_alert_success(
      "diamondback backend ready (Python {v$python}, NumPy {v$numpy}, SciPy {v$scipy})."
    )
  }
  invisible(mod)
}

#' Refuse to run an R half against a Python half it does not match
#'
#' The two are one program in two files, and they can drift apart. By far the
#' commonest way is a stale R session: the package is reinstalled while a session
#' still holds the old namespace, so the old R code calls the new Python module.
#' Without this check the first symptom is something like
#' `KeyError: 'edge_domain_len'` from deep inside a reduction, which is a
#' genuinely terrible way to be told to restart R.
#' @noRd
db_check_interface <- function(v, path) {
  got <- v$interface_version %||% "1"
  if (identical(as.character(got), PY_INTERFACE_VERSION)) {
    return(invisible(TRUE))
  }
  cli::cli_abort(c(
    "diamondback's R code and its Python module are out of step.",
    "x" = "The R code expects interface version {.val {PY_INTERFACE_VERSION}}; \\
           the Python module reports {.val {got}}.",
    "i" = "This usually means the package was reinstalled or updated while this R \\
           session had it loaded. R does not replace a namespace that is already \\
           in memory, so the R half is the old one and the Python half is the new one.",
    "*" = "Restart R and try again. In RStudio: {.emph Session > Restart R}, \\
           or {.kbd Ctrl/Cmd+Shift+F10}.",
    "*" = "If it persists after a restart, reinstall diamondback.",
    " " = "Python module: {.path {path}}"
  ), class = "diamondback_interface_mismatch", call = NULL)
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

#' A short, de-duplicated list of the environment names we can see
#' @noRd
db_known_env_names <- function(cands = NULL) {
  if (is.null(cands)) cands <- tryCatch(db_python_candidates(), error = function(e) list())
  if (!length(cands)) return(character())
  unique(vapply(cands, `[[`, character(1), "name"))
}

db_env_name_error <- function(name) {
  known <- db_known_env_names()
  bullets <- c(
    "No Python environment named {.val {name}} was found.",
    "x" = "{.code options(diamondback.python = {.val {name}})} matched neither a \\
           file path nor an environment name."
  )
  if (length(known)) {
    bullets <- c(
      bullets,
      "i" = "Environments diamondback can see: {.val {utils::head(known, 12)}}."
    )
  }
  bullets <- c(
    bullets,
    "i" = "You can also give a full path to a Python interpreter, or run \\
           {.run diamondback_check()} to see what was found."
  )
  cli::cli_abort(bullets, class = "diamondback_python_missing", call = NULL)
}

db_python_missing_error <- function() {
  cands <- tryCatch(db_python_candidates(), error = function(e) list())

  bullets <- c(
    "diamondback needs Python with NumPy and SciPy, and could not find one.",
    "x" = "No Python environment that diamondback can see has both installed."
  )

  # Leaving a user to hunt for a path when we have just enumerated every
  # interpreter on the machine would be unkind. Name them, so the fix is a
  # copy-paste rather than a search.
  if (length(cands)) {
    nm <- vapply(cands, `[[`, character(1), "name")
    py <- vapply(cands, `[[`, character(1), "python")
    shown <- utils::head(sprintf("%s  (%s)", nm, py), 8)
    bullets <- c(
      bullets,
      "i" = "Environments found, none of which has both NumPy and SciPy:",
      stats::setNames(shown, rep(" ", length(shown)))
    )
  }

  bullets <- c(
    bullets,
    "i" = "diamondback will not download anything on its own. Pick one:",
    "*" = "Install them into an environment you already have, e.g. \\
           {.code conda install -n <env> numpy scipy} or {.code pip install numpy scipy}.",
    "*" = "Point diamondback at an environment, by name or by path: \\
           {.code options(diamondback.python = \"myenv\")}.",
    "*" = "{.run diamondback_install_python()} -- set up a private environment \\
           for diamondback (downloads ~200 MB; asks first)."
  )
  cli::cli_abort(bullets, class = "diamondback_python_missing", call = NULL)
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
#'   `scipy`, `mode`, `interface_version`, `memory_gb`, `test_passed` and `ok`,
#'   invisibly.
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
    interface_version = NA_character_,
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
    out$interface_version <- v$interface_version %||% NA_character_
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
