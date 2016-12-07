#' @export
predict.cpls <- function(object, newdata, ncomp = 1:object$ncomp, comps,
                         type = c("response", "scores"), na.action = na.pass, ...) {
  x <- model.frame(object$formula, data = newdata)[, -1, drop=FALSE]
  pred <- t(object$alpha) %*% t(object$w) %*% t(x) + object$center
  as.vector(pred)
}
