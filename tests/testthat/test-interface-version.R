# The R code and inst/python/diamondback_core.py are one program in two files.
# They can drift apart -- most easily when a package is reinstalled while an R
# session still holds the old namespace, so the old R half calls the new Python
# half. The first symptom of that was "KeyError: 'edge_domain_len'" from inside
# a reduction, which is a terrible way to be told to restart R.

test_that("the R and Python interface versions agree", {
  skip_if_no_python()
  # The check that would have caught the drift at its source. If this fails,
  # PY_INTERFACE_VERSION and INTERFACE_VERSION have been bumped out of step.
  v <- py_r(diamondback:::db_py()$versions())
  expect_equal(v$interface_version, PY_INTERFACE_VERSION)
})

test_that("the declared interface version matches the shipped Python file", {
  # Reads the source rather than the loaded module, so a stale import cannot
  # make this pass.
  f <- system.file("python", "diamondback_core.py", package = "diamondback")
  if (!nzchar(f)) f <- "../../inst/python/diamondback_core.py"
  skip_if_not(file.exists(f), "Python source not found")

  txt <- readLines(f, warn = FALSE)
  line <- grep('^INTERFACE_VERSION\\s*=', txt, value = TRUE)
  expect_length(line, 1L)
  declared <- gsub('.*"(.*)".*', "\\1", line)
  expect_equal(declared, PY_INTERFACE_VERSION)
})

test_that("a mismatched Python module is refused with an actionable message", {
  # The scenario: R half expects one interface, Python half reports another.
  err <- tryCatch(
    db_check_interface(list(interface_version = "999"), "/some/path"),
    error = function(e) e
  )
  expect_s3_class(err, "diamondback_interface_mismatch")

  msg <- conditionMessage(err)
  expect_match(msg, "out of step")
  expect_match(msg, "999")                 # what was found
  expect_match(msg, PY_INTERFACE_VERSION)  # what was expected
  # The message must say what to actually do, since the cause is invisible.
  expect_match(msg, "Restart R")
})

test_that("an older Python module with no interface version is caught", {
  # versions() gained interface_version in interface 2; anything without it is
  # older than that and must not be run against this R code.
  expect_error(db_check_interface(list(), "/p"),
               class = "diamondback_interface_mismatch")
  expect_error(db_check_interface(list(python = "3.12"), "/p"), "out of step")
})

test_that("a matching interface passes silently", {
  expect_true(db_check_interface(list(interface_version = PY_INTERFACE_VERSION), "/p"))
  expect_silent(db_check_interface(list(interface_version = PY_INTERFACE_VERSION), "/p"))
})

test_that("diamondback_check reports the interface version", {
  skip_if_no_python()
  chk <- diamondback_check(verbose = FALSE)
  expect_equal(chk$interface_version, PY_INTERFACE_VERSION)
  expect_true(chk$ok)
})

test_that("every key the R code asks of the backend actually exists", {
  skip_if_no_python()
  # A direct guard against the class of bug that started this: an R call site
  # asking for a dict key the Python module no longer returns. Exercises the
  # real kernels and names any key that has gone missing.
  py <- diamondback:::db_py()
  code <- py$code_block(as.numeric(c(1, 1, 0, 0)), 2L, 2L, mask = NULL,
                        class_values = NULL, n_classes = 1L)
  lab <- py_get(py$label_array(code, 0L, 8L, FALSE), "labels")

  expected <- list(
    patch_stats = c("count", "sum_row", "sum_col"),
    patch_bboxes = c("row_min", "row_max", "col_min", "col_max"),
    edge_lengths = c("edge_valid_len", "edge_missing_len", "edge_outside_len",
                     "edge_ns_cells", "edge_ew_cells", "edge_valid_cells",
                     "edge_missing_cells", "edge_outside_cells"),
    core_counts = c("core_count", "core_mask"),
    overlap_counts = c("id1", "id2", "cells")
  )
  got <- list(
    patch_stats = py$patch_stats(lab, 1L),
    patch_bboxes = py$patch_bboxes(lab, 1L),
    edge_lengths = py$edge_lengths(lab, code, 1L, c(1, 1), c(1, 1, 1), FALSE),
    core_counts = py$core_counts(code, lab, 1L, 0L, 0.5,
                                 reticulate::tuple(1, 1), "all"),
    overlap_counts = py$overlap_counts(lab, lab, 1L, 1L)
  )

  for (fn in names(expected)) {
    keys <- names(py_r(got[[fn]]))
    missing <- setdiff(expected[[fn]], keys)
    expect_length(missing, 0L)
    if (length(missing)) {
      fail(sprintf("%s() no longer returns: %s", fn, paste(missing, collapse = ", ")))
    }
  }
})
