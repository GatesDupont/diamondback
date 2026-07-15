# Finding a Python should not require the user to know where one lives. These
# tests cover the resolution of environment *names*, and the auto-detection that
# means most users never set anything at all.

fake_env <- function(name, dir = withr::local_tempdir(.local_envir = parent.frame())) {
  # A directory shaped like an environment prefix, with an executable stub where
  # a real interpreter would be.
  bin <- if (.Platform$OS.type == "windows") "Scripts" else "bin"
  d <- file.path(dir, "envs", name, bin)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  exe <- file.path(d, if (.Platform$OS.type == "windows") "python.exe" else "python")
  file.create(exe)
  Sys.chmod(exe, "0755")
  list(root = dir, name = name, python = exe)
}

test_that("db_python_bin points at the interpreter inside a prefix", {
  p <- db_python_bin("/somewhere/env")
  expect_true(any(grepl("python", basename(p))))
  if (.Platform$OS.type != "windows") {
    expect_equal(p, file.path("/somewhere/env", "bin", "python"))
  }
})

test_that("conda roots are found from CONDA_EXE and CONDA_PREFIX", {
  d <- withr::local_tempdir()
  dir.create(file.path(d, "bin"), recursive = TRUE)
  file.create(file.path(d, "bin", "conda"))

  withr::local_envvar(CONDA_EXE = file.path(d, "bin", "conda"), CONDA_PREFIX = "")
  expect_true(normalizePath(d) %in% normalizePath(db_conda_roots(), mustWork = FALSE))

  # An *active* environment implies the root two levels up.
  envd <- file.path(d, "envs", "geo")
  dir.create(envd, recursive = TRUE)
  withr::local_envvar(CONDA_EXE = "", CONDA_PREFIX = envd)
  roots <- normalizePath(db_conda_roots(), mustWork = FALSE)
  expect_true(normalizePath(d) %in% roots)      # the root above the env
  expect_true(normalizePath(envd) %in% roots)   # and the env itself
})

test_that("an environment is found by name, without any path", {
  e <- fake_env("geo")
  withr::local_envvar(CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  local_mocked_bindings(db_conda_roots = function() e$root)

  expect_equal(normalizePath(db_resolve_env_name("geo")), normalizePath(e$python))
  expect_null(db_resolve_env_name("no_such_env"))
})

test_that("environment names are listed for the error message", {
  e <- fake_env("geo")
  withr::local_envvar(CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  local_mocked_bindings(db_conda_roots = function() e$root)
  expect_true("geo" %in% db_known_env_names())
})

test_that("naming an environment resolves it and sets RETICULATE_PYTHON", {
  local_fresh <- function() testthat::local_mocked_bindings(
    db_py_initialized = function() FALSE, .env = parent.frame(2))
  local_fresh()
  e <- fake_env("geo")
  withr::local_envvar(RETICULATE_PYTHON = "", DIAMONDBACK_PYTHON = "",
                      CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  withr::local_options(diamondback.python = "geo")
  local_mocked_bindings(
    db_conda_roots = function() e$root,
    db_has_deps = function(python) TRUE          # stub has no real interpreter
  )

  expect_equal(db_python_mode(), "user")
  expect_equal(normalizePath(Sys.getenv("RETICULATE_PYTHON")), normalizePath(e$python))
})

test_that("a name that matches nothing lists the names that do exist", {
  testthat::local_mocked_bindings(db_py_initialized = function() FALSE)
  e <- fake_env("geo")
  withr::local_envvar(RETICULATE_PYTHON = "", DIAMONDBACK_PYTHON = "",
                      CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  withr::local_options(diamondback.python = "nope")
  local_mocked_bindings(db_conda_roots = function() e$root)

  err <- tryCatch(db_python_mode(), error = function(e) e)
  expect_s3_class(err, "diamondback_python_missing")
  msg <- conditionMessage(err)
  expect_match(msg, "No Python environment named")
  expect_match(msg, "geo")           # tells you what you could have said
})

test_that("a path-shaped value that does not exist still reads as a bad path", {
  # Distinct from a bad *name*: the message should not offer a list of names.
  testthat::local_mocked_bindings(db_py_initialized = function() FALSE)
  withr::local_envvar(RETICULATE_PYTHON = "", DIAMONDBACK_PYTHON = "")
  withr::local_options(diamondback.python = "/no/such/python")
  expect_error(db_python_mode(), "does not exist")
})

test_that("a named environment without NumPy or SciPy is refused clearly", {
  testthat::local_mocked_bindings(db_py_initialized = function() FALSE)
  e <- fake_env("bare")
  withr::local_envvar(RETICULATE_PYTHON = "", DIAMONDBACK_PYTHON = "",
                      CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  withr::local_options(diamondback.python = "bare")
  local_mocked_bindings(
    db_conda_roots = function() e$root,
    db_has_deps = function(python) FALSE
  )

  err <- tryCatch(db_python_mode(), error = function(e) e)
  expect_s3_class(err, "diamondback_python_missing")
  expect_match(conditionMessage(err), "missing NumPy or SciPy")
})

test_that("auto-detection picks the first environment that has the dependencies", {
  e1 <- fake_env("aaa"); e2 <- fake_env("zzz", dir = e1$root)
  withr::local_envvar(CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  local_mocked_bindings(
    db_conda_roots = function() e1$root,
    # Only the second one qualifies, so the first must be passed over rather
    # than picked because it came first.
    db_has_deps = function(python) grepl("zzz", python)
  )
  expect_equal(normalizePath(db_probe_python()), normalizePath(e2$python))
})

test_that("auto-detection is deterministic across calls", {
  e1 <- fake_env("aaa"); fake_env("mmm", dir = e1$root); fake_env("zzz", dir = e1$root)
  withr::local_envvar(CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  local_mocked_bindings(
    db_conda_roots = function() e1$root,
    # Only the fake environments qualify. Without this the real python3 on the
    # PATH would win -- correctly, since PATH is a stronger signal than a
    # scanned environment, but that is a different test than this one.
    db_has_deps = function(python) grepl(e1$root, python, fixed = TRUE)
  )
  # Environments are sorted, so the same machine always gives the same answer.
  expect_equal(db_probe_python(), db_probe_python())
  expect_match(db_probe_python(), "aaa")
})

test_that("a Python on the PATH is preferred over a scanned environment", {
  # Ordering is a design choice, not an accident: what `python` already means
  # here is a stronger signal than an environment we went looking for.
  e <- fake_env("aaa")
  withr::local_envvar(CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  local_mocked_bindings(
    db_conda_roots = function() e$root,
    db_has_deps = function(python) TRUE
  )
  skip_if(!nzchar(Sys.which("python3")), "no python3 on PATH")
  expect_equal(normalizePath(db_probe_python()),
               normalizePath(unname(Sys.which("python3"))))
})

test_that("candidate discovery survives a broken conda_list()", {
  # reticulate::conda_list() returns a bogus row on some miniforge layouts, and
  # errors on machines with no conda at all. Neither may take discovery down.
  e <- fake_env("geo")
  withr::local_envvar(CONDA_EXE = "", CONDA_PREFIX = "", VIRTUAL_ENV = "", WORKON_HOME = "")
  local_mocked_bindings(db_conda_roots = function() e$root)

  testthat::local_mocked_bindings(
    conda_list = function(...) stop("no conda"), .package = "reticulate")
  expect_true("geo" %in% db_known_env_names())

  testthat::local_mocked_bindings(
    conda_list = function(...) data.frame(name = "envs", python = "/does/not/exist"),
    .package = "reticulate")
  names_found <- db_known_env_names()
  expect_true("geo" %in% names_found)
  expect_false("envs" %in% names_found)   # the bogus entry is dropped
})

test_that("db_has_deps accepts a real interpreter and rejects a fake one", {
  skip_if_no_python()
  cfg <- reticulate::py_config()
  expect_true(db_has_deps(cfg$python))
  expect_false(db_has_deps(file.path(tempdir(), "definitely_not_python")))
})
