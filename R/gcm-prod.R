gcm_prod <- function(matrix1, matrix2) {
  UseMethod("gcm_prod")
}

gcm_prod.default <- `%*%`

gcm_prod.GrammarCompressedMatrix <- function(gcmatrix, matrix) {
  gcmatrix$prod(matrix)
}

gcm_tprod <- function(matrix1, matrix2) {
  UseMethod("gcm_tprod")
}

gcm_tprod.default <- function(matrix1, matrix2) {
  t(matrix1) %*% matrix2
}

gcm_tprod.GrammarCompressedMatrix <- function(gcmatrix, matrix) {
  gcmatrix$tprod(matrix)
}
