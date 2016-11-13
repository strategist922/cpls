gcmprod <- function(matrix1, matrix2) {
  UseMethod("gcmprod")
}

gcmprod.default <- `%*%`

gcmprod.GrammarCompressedMatrix <- function(gcmatrix, matrix) {
  gcmatrix$prod(matrix)
}

gcmtprod <- function(matrix1, matrix2) {
  UseMethod("gcmtprod")
}

gcmtprod.default <- function(matrix1, matrix2) {
  t(matrix1) %*% matrix2
}

gcmtprod.GrammarCompressedMatrix <- function(gcmatrix, matrix) {
  gcmatrix$tprod(matrix)
}
