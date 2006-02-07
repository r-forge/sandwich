## bread.default <- function(x, ...)
## {
##   X <- estfunDeriv(x, ...)
##   rval <- apply(X, 1:2, sum)/dim(X)[3]
##   rval <- solve(rval)
##   rownames(rval) <- colnames(rval) <- dimnames(X)[[1]]
##   return(rval)
## }

## meat <- function(x, ...)
## {
##   UseMethod("meat")
## }

## meat.default -> see meat

estfunDeriv <- function(x, ...)
{
  UseMethod("estfunDeriv")
}

estfunDeriv.lm <- function(x, ...)
{
  ## extract X
  xmat <- if(is.matrix(x$x)) x$x else model.matrix(x)
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

## check: glm & survreg
## still missing: coxph

estfunDeriv.glm <- function(x, ...)
{
  ## extract X from IWLS fit
  xmat <- qr.X(x$qr)
  attr(xmat, "assign") <- NULL
  
  ## compute cross-products for each observation i
  k <- NCOL(xmat)
  n <- NROW(xmat)
  rval <- sapply(1:n, function(i) crossprod(xmat[i,,drop = FALSE]))
  dim(rval) <- c(k, k, n)
  dimnames(rval) <- list(colnames(xmat), colnames(xmat), rownames(xmat))
    
  return(rval)
}

estfunDeriv.survreg <- function(x, ...)
{
  stopifnot(require(survival))
  xmat <- if(is.matrix(x$x)) x$x else model.matrix(x)
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

