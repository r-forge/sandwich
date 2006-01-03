sandwich <- function(x, bread. = bread, meat. = meat, ...)
{
  if(is.function(bread.)) bread. <- bread.(x)
  if(is.function(meat.)) meat. <- meat.(x, ...)
  n <- NROW(estfun(x))
  return(1/n * (bread. %*% meat. %*% bread.))
}

bread <- function(x, ...)
{
  UseMethod("bread")
}

## bread.default <- function(x, ...)
## {
##   X <- estfunDeriv(x, ...)
##   rval <- apply(X, 1:2, sum)/dim(X)[3]
##   rval <- solve(rval)
##   rownames(rval) <- colnames(rval) <- dimnames(X)[[1]]
##   return(rval)
## }

bread.lm <- bread.nls <- function(x, ...)
{
  sx <- summary(x)
  return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])))
}

bread.survreg <- function(x, ...)
{
  stopifnot(require("survival"))

  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  wts <- if(!is.null(x$weights)) x$weights else 1
  xmat <- as.matrix(wts * xmat)
  attr(xmat, "assign") <- NULL

  n <- NROW(xmat)
  res <- residuals(x, type = "matrix")

  if(NROW(x$var) > length(coef(x))) {
    rval <- crossprod(xmat, res[,"ddg"] * xmat)/n
    xmean <- colMeans(res[,"dsg"] * xmat)
    rval <- cbind(rbind(rval, xmean), c(xmean, mean(res[,"dds"])))
    rownames(rval) <- colnames(rval) <- c(colnames(xmat), "Log(scale)")
  } else {
    rval <- crossprod(xmat, res[,"ddg"] * xmat)/n
    rownames(rval) <- colnames(rval) <- colnames(xmat)
  }

  rval <- solve(rval)
  return(rval)
}

## meat <- function(x, ...)
## {
##   UseMethod("meat")
## }

## meat.default <- function(x, adjust = FALSE, ...)
meat <- function(x, adjust = FALSE, ...)
{
  psi <- estfun(x, ...)
  k <- NCOL(psi)
  n <- NROW(psi)
  rval <- crossprod(as.matrix(psi))/n
  if(adjust) rval <- n/(n-k) * rval
  rownames(rval) <- colnames(rval) <- colnames(psi)
  return(rval)
}
