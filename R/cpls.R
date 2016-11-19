#' Compression-based Partial Least Squares Regression (cPLS)
#'
#' @param formula formula
#' @param data data
#' @param ncomp ncomp
#' @param verbose verbose
#'
#' @export
cpls <- function(formula, data = NULL, headers = NULL, ncomp = 10,
                 verbose = TRUE) {
  UseMethod("cpls", object = data)
}

#' @importFrom stats model.frame
#' @export
cpls.default <- function(formula, data, headers, ncomp, verbose) {
  if (is.null(data)) {
    data <- parent.frame()
  }
  df <- model.frame(formula = formula, data = data)
  y <- df[, 1]
  x <- as.matrix(df[, -1])
  cpls.matrix(y, x, headers, ncomp, verbose)
}

cpls.character <- function(formula, data, headers, ncomp) {

}

cpls.matrix <- function(y, x, headers, ncomp, verbose) {
  y <- scale(y, scale = FALSE)
  center <- attr(y, "scaled:center")
  y <- matrix(y)
  gcm <- GrammarCompressedMatrix$new(x, verbose = verbose)
  cpls.GrammarCompressedMatrix(y, gcm, ncomp, center, verbose)
}

cpls.GrammarCompressedMatrix <- function(y, gcm, ncomp, center, verbose) {
  d <- gcm$ncol()
  n <- gcm$nrow()
  w <- matrix(0, nrow = d, ncol = ncomp)
  t <- matrix(0, nrow = n, ncol = d)
  r <- y
  for (i in seq_len(ncomp)) {
    w[,i] <- gcm_tprod(gcm, r)
    if (i == 1) {
      t[,1] <- gcm_prod(gcm, w[,1])
    } else {
      t[,i] <- (diag(n) - t[,seq_len(i-1)] %*% t(t[,seq_len(i-1)])) %*% gcm_prod(gcm, w[,1])
    }
    t[,i] <- t[,i] / sqrt(sum(t[,i]^2))
    r <- r - as.numeric(t(y) %*% t[,i]) * t[,i,drop=FALSE]
  }
  # alpha <- solve(t(w) %*% t(x) %*% x %*% w, t(w) %*% t(x) %*% y)
  alpha <- solve(t(w) %*% gcm_tprod(gcm, gcm_prod(gcm, w)), t(w) %*% gcm_tprod(gcm, y))
  list(alpha=alpha, w=w, x=gcm, y=y, center=center)
}
