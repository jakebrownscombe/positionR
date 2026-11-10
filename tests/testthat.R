# testthat is in Suggests, and the hard-dependencies-only CI job checks without it.
if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(positionR)
  test_check("positionR")
}
