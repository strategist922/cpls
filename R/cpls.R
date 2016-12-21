#' Compression-based Partial Least Squares Regression (cPLS)
#'
#' @param formula formula
#' @param data data
#' @param ncomp ncomp
#' @param re_pair_method re_pair_method
#' @param verbose verbose
#'
#' @importFrom stats model.frame
#' @export
cpls <- function(formula, data = NULL, ncomp = 10,
                 re_pair_method = c("freq", "naive", "lossless", "lossy"),
                 verbose = TRUE) {
  re_pair_method <- match.arg(re_pair_method)
  df <- model.frame(formula = formula, data = data)
  y <- scale(df[,1], scale = FALSE)
  center <- attr(y, "scaled:center")
  y <- matrix(y)
  x <- as.matrix(df[, -1, drop=FALSE])
  d <- ncol(x)
  n <- nrow(x)
  gcm <- GrammarCompressedMatrix$new(x, re_pair_method = re_pair_method, verbose = verbose)

  if (verbose) message("################## cPLS #################################")
  if (verbose) progress_bar <- txtProgressBar(0, ncomp+1, style=3)
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
    if (verbose) setTxtProgressBar(progress_bar, i)
  }
  alpha <- solve(t(w) %*% gcm_tprod(gcm, gcm_prod(gcm, w)), t(w) %*% gcm_tprod(gcm, y))
  if (verbose) setTxtProgressBar(progress_bar, ncomp+1)
  result <- list(formula=formula, ncomp=ncomp, alpha=alpha, w=w, center=center)
  class(result) <- "cpls"
  if (verbose) message("################## Finished #############################")
  result
}
