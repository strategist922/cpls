context("gcm_prod")

test_that("prod for GCM", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm_prod(gcm, mat), x %*% mat)
})

test_that("tprod for GCM", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm_tprod(gcm, mat), t(x) %*% mat)
})

test_that("prod for matrix", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm_prod(x, mat), x %*% mat)
})

test_that("tprod for matrix", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm_tprod(x, mat), t(x) %*% mat)
})
