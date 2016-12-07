context("GrammarCompressedMatrix")

test_that("Create Object", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  obj <- GrammarCompressedMatrix$new(x, re_pair_method ="freq", verbose = FALSE)
  expect_identical(class(obj), c("GrammarCompressedMatrix", "R6"))
})

test_that("Uncompress Row", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method ="freq", verbose = FALSE)
  for (i in seq_len(nrow(x))) {
    expect_identical(gcm$access_row(i), x[i, ])
  }
})

test_that("Uncompress Column", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method ="freq", verbose = FALSE)
  for (i in seq_len(ncol(x))) {
    expect_identical(gcm$access_column(i), x[, i])
  }
})

test_that("prod", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method ="freq", verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm$prod(mat), x %*% mat)
})

test_that("tprod", {
  x <- matrix(sample(0:1, size = 100, replace = TRUE), nrow=10)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method ="freq", verbose = FALSE)
  mat <- matrix(1:10)
  expect_identical(gcm$tprod(mat), t(x) %*% mat)
})
