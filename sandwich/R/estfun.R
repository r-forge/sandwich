estfun <- function(x, ...)
{
  UseMethod("estfun")
}

estfun.lm <- function(x, ...)
{
  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  wts <- if(!is.null(x$weights)) x$weights else 1
  res <- residuals(x)
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(is.zoo(res)) rval <- zoo(rval, time(res))
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  return(rval)
}

estfun.glm <- function(x, ...)
{
  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  rval <- as.vector(residuals(x, "working")) * x$weights * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  res <- residuals(x, type = "pearson")
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(is.zoo(res)) rval <- zoo(rval, time(res))
  return(rval)
}

estfun.rlm <- function(x, ...)
{
  if (is.matrix(x$x)) 
      xmat <- x$x
  else {
      mf <- model.frame(x)
      xmat <- model.matrix(terms(x), mf)
  }
  if (!is.null(x$weights)) 
      wts <- x$weights
  else wts <- 1
  res <- residuals(x)
  psi <- function(z) x$psi(z) * z
  rval <- as.vector(psi(res/x$s)) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
  if(is.zoo(res)) rval <- zoo(rval, time(res))
  return(rval)
}

######################################################
## new ###############################################
######################################################

estfunDeriv <- function(x, ...)
{
  UseMethod("estfunDeriv")
}

estfunDeriv.lm <- function(x, ...)
{
  ## extract X
  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  wts <- if(!is.null(x$weights)) x$weights else 1
  xmat <- as.matrix(wts * xmat)
  attr(xmat, "assign") <- NULL
  
  ## compute cross-products for each observation i
  k <- NCOL(xmat)
  n <- NROW(xmat)
  rval <- sapply(1:n, function(i) crossprod(xmat[i,,drop = FALSE]))
  dim(rval) <- c(k, k, n)
  dimnames(rval) <- list(colnames(xmat), colnames(xmat), rownames(xmat))
    
  return(rval)
}

bread <- function(x, ...)
{
  UseMethod("bread")
}

bread.default <- function(x, ...)
{
  X <- estfunDeriv(x, ...)
  rval <- apply(X, 1:2, sum)/dim(X)[3]
  rval <- solve(rval)
  rownames(rval) <- colnames(rval) <- dimnames(X)[[1]]
  return(rval)
}

bread.lm <- function(x, ...)
{
  sx <- summary(x)
  return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])))
}

meat <- function(x, ...)
{
  UseMethod("meat")
}

meat.default <- function(x, adjust = FALSE, ...)
{
  psi <- estfun(x, ...)
  k <- NCOL(psi)
  n <- NROW(psi)
  rval <- crossprod(as.matrix(psi))/(n^2)
  if(adjust) rval <- n/(n-k) * rval
  rownames(rval) <- colnames(rval) <- colnames(psi)
  return(rval)
}

sandwich <- function(x, breadfun = bread, meatfun = meat, ...)
{
  rval <- breadfun(x)
  rval <- rval %*% meatfun(x, ...) %*% rval
  return(rval)  
}


## survival functions

estfun.coxph <- function(x, ...)
{
  stopifnot(require(survival))
  residuals(x, type = "score", ...)
}

estfun.survreg <- function(x, ...)
{
  stopifnot(require(survival))
  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  wts <- if(!is.null(x$weights)) x$weights else 1
  res <- residuals(x, type = "matrix")
  rval <- as.vector(res[,"dg"]) * wts * xmat * 2 #FIXME# like this?
  if(NROW(x$var) > length(coefficients(x))) {
    rval <- cbind(rval, res[,"ds"])
    colnames(rval)[NCOL(rval)] <- "Log(scale)"
  }
  attr(rval, "assign") <- NULL
  
  return(rval)
}

estfunDeriv.survreg <- function(x, ...)
{
  stopifnot(require(survival))
  if (is.matrix(x$x))
    xmat <- x$x
  else {
    mf <- model.frame(x)
    xmat <- model.matrix(terms(x), mf)    
  }
  wts <- if(!is.null(x$weights)) x$weights else 1
  xmat <- as.matrix(wts * xmat)
  attr(xmat, "assign") <- NULL

  k <- NCOL(xmat)
  n <- NROW(xmat)
  res <- residuals(x, type = "matrix")

  if(NROW(x$var) > length(coefficients(x))) {
    funi <- function(i) {
      xi <- xmat[i,,drop = FALSE]
      resi <- res[i,]
      rvali <- crossprod(resi["ddg"] * xi)
      rvali <- rbind(rvali, resi["dsg"] * xi)
      rvali <- cbind(rvali, c(resi["dsg"] * as.vector(xi), resi["dds"]))
      return(rvali)      
    }
    rval <- sapply(1:n, funi)
    dim(rval) <- c(k+1, k+1, n)
    dimnames(rval) <- list(c(colnames(xmat), "Log(scale)"), c(colnames(xmat), "Log(scale)"), rownames(xmat))
  } else {
    res <- res[,"ddg"]
    rval <- sapply(1:n, function(i) crossprod(res[i] * xmat[i,,drop = FALSE]))
    dim(rval) <- c(k, k, n)
    dimnames(rval) <- list(colnames(xmat), colnames(xmat), rownames(xmat))
  }

  return(rval)
}

