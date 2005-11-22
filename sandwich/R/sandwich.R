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
