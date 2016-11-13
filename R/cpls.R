#' @export
cpls <- function(formula, data = NULL, ncomp = 10) {
  df <- model.frame(formula = formula, data = data)
  y <- scale(df[,1], scale = FALSE)
  center <- attr(y, "scaled:center")
  y <- matrix(y)
  x <- as.matrix(df[, -1])
  d <- ncol(x)
  n <- nrow(x)
  gcm <- GrammarCompressedMatrix$new(x)

  w <- matrix(0, nrow = d, ncol = ncomp)
  t <- matrix(0, nrow = n, ncol = d)
  r <- y
  for (i in seq_len(ncomp)) {
    # w[,i] <- t(x) %*% r
    w[,i] <- gcmtprod(gcm, r)
    if (i == 1) {
      # t[,1] <- x %*% w[,1]
      t[,1] <- gcmprod(gcm, w[,1])
    } else {
      # t[,i] <- (diag(n) - t[,seq_len(i-1)] %*% t(t[,seq_len(i-1)])) %*% x %*% w[,i]
      t[,i] <- (diag(n) - t[,seq_len(i-1)] %*% t(t[,seq_len(i-1)])) %*% gcmprod(gcm, w[,1])
    }
    t[,i] <- t[,i] / sqrt(sum(t[,i]^2))
    r <- r - as.numeric(t(y) %*% t[,i]) * t[,i,drop=FALSE]
  }
  alpha <- solve(t(w) %*% t(x) %*% x %*% w, t(w) %*% t(x) %*% y)
  list(alpha=alpha, w=w, x=x, y=y, center=center)
}
