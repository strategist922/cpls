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
cpls.default <- function(formula, data, headers, ncomp, verbose) {
  if (is.null(data)) {
    data <- parent.frame()
  }
  df <- model.frame(formula = formula, data = data)

}

cpls.character <- function(formula, data, headers, ncomp) {

}

cpls.matrix <- function(y, x, headers, ncomp, verbose) {
  y <- scale(df[,1], scale = FALSE)
  center <- attr(y, "scaled:center")
  y <- matrix(y)
  x <- as.matrix(df[, -1])
  d <- ncol(x)
  n <- nrow(x)
  gcm <- GrammarCompressedMatrix$new(x, verbose = verbose)

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
  alpha <- solve(t(w) %*% t(x) %*% x %*% w, t(w) %*% t(x) %*% y)
  list(alpha=alpha, w=w, x=x, y=y, center=center)
}
