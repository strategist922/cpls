context("gcm_prod")

test_that("prod for GCM", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method ="plain", verbose = FALSE)
  mat <- matrix(1:10)
  # print(x[c(2,7,9,10),])
  # print(gcm$access_row(2))
  # print(gcm$access_row(7))
  # print(gcm$access_row(9))
  # print(gcm$access_row(10))
  # print(gcm$.__enclos_env__$private$.compressed_matrix[c(2,7,9,10)])
  # print(gcm$.__enclos_env__$private$.grammer_trees)
  expect_identical(gcm_prod(gcm, mat), x %*% mat)
})

test_that("tprod for GCM", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method ="freq", verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm_tprod(gcm, mat), t(x) %*% mat)
})

test_that("prod for matrix", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method ="freq", verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm_prod(x, mat), x %*% mat)
})

test_that("tprod for matrix", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method ="freq", verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm_tprod(x, mat), t(x) %*% mat)
})
