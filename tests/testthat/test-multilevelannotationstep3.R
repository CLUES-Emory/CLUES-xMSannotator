test_that("HMDB data load does not pollute .GlobalEnv", {
  # Snapshot .GlobalEnv names before calling the internal loader
  before <- ls(envir = .GlobalEnv)

  # Load hmdbAllinf the same way the fixed multilevelannotationstep3 does:
  # envir = environment() keeps it local, not in .GlobalEnv
  data(hmdbAllinf, envir = environment())

  after <- ls(envir = .GlobalEnv)

  # No new objects should appear in .GlobalEnv
  expect_equal(before, after)
  # The object should exist locally
  expect_true(exists("hmdbAllinf", inherits = FALSE))
})
