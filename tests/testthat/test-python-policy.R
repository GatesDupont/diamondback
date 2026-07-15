# The promise is that running an analysis never downloads or installs anything.
# These tests pin the decision logic that enforces it. They deliberately test
# db_python_mode(), which decides *without* touching an interpreter -- that
# separation is the whole mechanism.

# Once any test initialises Python, db_python_mode() short-circuits to
# "existing" and the branches worth testing never run. Pretending Python is not
# yet up is what makes the decision logic testable at all.
local_fresh_session <- function(env = parent.frame()) {
  testthat::local_mocked_bindings(db_py_initialized = function() FALSE, .env = env)
}

test_that("an explicit RETICULATE_PYTHON is honoured and never installed over", {
  local_fresh_session()
  withr::local_envvar(RETICULATE_PYTHON = "/some/python", DIAMONDBACK_PYTHON = "")
  withr::local_options(diamondback.python = NULL)
  expect_equal(db_python_mode(), "user")
})

test_that("a configured path that does not exist is an error, not a download", {
  local_fresh_session()
  withr::local_envvar(RETICULATE_PYTHON = "", DIAMONDBACK_PYTHON = "")
  withr::local_options(diamondback.python = "/no/such/python")
  expect_error(db_python_mode(), "does not exist")
})

test_that("the managed environment is used only after explicit consent", {
  local_fresh_session()
  withr::local_envvar(RETICULATE_PYTHON = "", DIAMONDBACK_PYTHON = "")
  withr::local_options(diamondback.python = "managed")
  expect_equal(db_python_mode(), "managed")
})

test_that("consent recorded on disk enables the managed environment", {
  local_fresh_session()
  withr::local_envvar(RETICULATE_PYTHON = "", DIAMONDBACK_PYTHON = "")
  withr::local_options(diamondback.python = NULL)

  # Without consent and without a suitable system Python, the answer must be
  # "none" -- i.e. stop and ask -- rather than "managed".
  local_mocked_bindings(
    db_probe_python = function() NULL,
    db_consent_file = function() file.path(tempdir(), "no_such_consent_file")
  )
  expect_equal(db_python_mode(), "none")

  # With consent present, the managed environment becomes available.
  f <- withr::local_tempfile()
  writeLines("consented", f)
  local_mocked_bindings(
    db_probe_python = function() NULL,
    db_consent_file = function() f
  )
  expect_equal(db_python_mode(), "managed")
})

test_that("a system Python with NumPy and SciPy is used without any download", {
  local_fresh_session()
  withr::local_envvar(RETICULATE_PYTHON = "", DIAMONDBACK_PYTHON = "")
  withr::local_options(diamondback.python = NULL)

  local_mocked_bindings(
    db_probe_python = function() "/pretend/python3",
    db_consent_file = function() file.path(tempdir(), "no_such_consent_file")
  )
  expect_equal(db_python_mode(), "system")
  expect_equal(Sys.getenv("RETICULATE_PYTHON"), "/pretend/python3")
})

test_that("with no Python and no consent, the error explains rather than installs", {
  err <- tryCatch(db_python_missing_error(), error = function(e) e)
  expect_s3_class(err, "diamondback_python_missing")

  msg <- conditionMessage(err)
  expect_match(msg, "will not download anything on its own")
  expect_match(msg, "diamondback_install_python")
  expect_match(msg, "pip install numpy scipy")
})

test_that("the probe only accepts a Python that has both NumPy and SciPy", {
  # The real probe shells out; this checks it does not accept a bare interpreter.
  # Our own system python3 has NumPy but (on this machine) not SciPy.
  p <- db_probe_python()
  if (!is.null(p)) {
    out <- suppressWarnings(system2(p, c("-c", shQuote("import numpy, scipy")),
                                    stdout = TRUE, stderr = TRUE))
    expect_true(is.null(attr(out, "status")) || identical(attr(out, "status"), 0L))
  } else {
    succeed("no system Python with NumPy and SciPy, which is a valid outcome")
  }
})

test_that("diamondback_install_python declines cleanly when refused", {
  # menu() returning 2 is "No". Nothing should be downloaded or recorded.
  f <- file.path(tempdir(), "consent_should_not_appear")
  unlink(f)
  local_mocked_bindings(
    db_consent_file = function() f
  )
  testthat::local_mocked_bindings(menu = function(...) 2L, .package = "utils")

  expect_message(res <- diamondback_install_python(prompt = TRUE), "Cancelled")
  expect_false(res)
  expect_false(file.exists(f))
})

test_that("diamondback_check reports which Python source was used", {
  skip_if_no_python()
  chk <- diamondback_check(verbose = FALSE)
  expect_true(chk$ok)
  expect_true(chk$mode %in% c("managed", "system", "user", "existing"))
})
