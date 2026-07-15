# The guard exists to prevent a surprise OOM, not to prevent work. The first
# design refused whenever the estimate exceeded a fraction of *free* memory,
# which on macOS meant reading 1.7 GB free on an 8 GB machine and blocking a
# 4 GB job that runs there routinely. These tests pin the corrected division of
# labour: physical RAM may refuse, free memory may only warn.

test_that("physical memory is detected and exceeds free memory", {
  phys <- db_physical_memory()
  skip_if(is.na(phys), "physical memory not detectable on this platform")
  expect_gt(phys, 0)
  expect_gt(phys, 1024^3)          # any machine running this has > 1 GB

  free <- db_available_memory()
  if (!is.na(free)) {
    # The whole point: free is a snapshot and is smaller, often far smaller.
    expect_lte(free, phys)
  }
})

test_that("a job that cannot fit in physical RAM is refused", {
  expect_error(
    db_check_memory(1e12, "label", memory_limit = 8 * 1024^3, quiet = TRUE),
    class = "diamondback_memory_error"
  )
  err <- tryCatch(db_check_memory(1e12, "label", memory_limit = 8 * 1024^3, quiet = TRUE),
                  error = conditionMessage)
  expect_match(err, "more than this machine has")
  expect_match(err, "structural limit")     # says why it is a refusal, not a guess
})

test_that("a job that fits physical RAM is allowed even when free memory is low", {
  # This is the regression. 4 GB needed, 8 GB physical, only 1 GB free: the old
  # guard refused; the new one runs and warns.
  local_mocked_bindings(db_available_memory = function() 1 * 1024^3)
  expect_no_error(
    db_check_memory(724461424, "label", memory_limit = 8 * 1024^3, quiet = TRUE)
  )
})

test_that("low free memory warns rather than refusing", {
  local_mocked_bindings(
    db_available_memory = function() 1 * 1024^3,
    db_physical_memory = function() 8 * 1024^3
  )
  expect_warning(db_check_memory(724461424, "label", quiet = FALSE), NA)  # not a warning()
  msgs <- capture.output(db_check_memory(724461424, "label", quiet = FALSE), type = "message")
  txt <- paste(msgs, collapse = " ")
  expect_match(txt, "free right now")
  expect_match(txt, "paging")          # tells the user what to expect
})

test_that("a comfortable job says nothing alarming", {
  local_mocked_bindings(
    db_available_memory = function() 32 * 1024^3,
    db_physical_memory = function() 64 * 1024^3
  )
  msgs <- capture.output(db_check_memory(1e6, "label", quiet = FALSE), type = "message")
  txt <- paste(msgs, collapse = " ")
  expect_match(txt, "Estimated peak memory")
  expect_no_match(txt, "free right now")
  expect_no_match(txt, "large share")
})

test_that("memory_limit replaces physical RAM, not free memory", {
  local_mocked_bindings(
    db_physical_memory = function() 64 * 1024^3,   # a big machine ...
    db_available_memory = function() 32 * 1024^3
  )
  # ... but the caller says the real budget is 1 GB, e.g. a container cgroup.
  expect_error(
    db_check_memory(1e9, "label", memory_limit = 1 * 1024^3, quiet = TRUE),
    class = "diamondback_memory_error"
  )
})

test_that("undetectable physical memory warns but does not block", {
  local_mocked_bindings(
    db_physical_memory = function() NA_real_,
    db_available_memory = function() NA_real_
  )
  expect_no_error(db_check_memory(1e6, "label", quiet = TRUE))
  expect_warning(db_check_memory(4e9, "label", quiet = TRUE), "Could not detect physical memory")
})

test_that("max_memory_frac scales the hard ceiling", {
  local_mocked_bindings(db_available_memory = function() 8 * 1024^3)
  # 5 GB needed against 8 GB physical: allowed at 0.9, refused at 0.5.
  expect_no_error(db_check_memory(5e9 / 6, "label", max_memory_frac = 0.9,
                                  memory_limit = 8 * 1024^3, quiet = TRUE))
  expect_error(db_check_memory(5e9 / 6 * 6, "label", max_memory_frac = 0.05,
                               memory_limit = 8 * 1024^3, quiet = TRUE),
               class = "diamondback_memory_error")
})

test_that("no hand-tuned constant is needed for a Central-Hardwoods-sized job", {
  skip_if(is.na(db_physical_memory()), "physical memory not detectable")
  skip_if(db_physical_memory() < 6 * 1024^3, "test assumes >= 6 GB of RAM")
  # 724M cells, labelling and core area, with the real detector and no override.
  expect_no_error(db_check_memory(724461424, "label", quiet = TRUE))
  expect_no_error(db_check_memory(724461424, "core", quiet = TRUE))
})
