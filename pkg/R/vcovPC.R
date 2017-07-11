vcovPC <- function(x, cluster = NULL, order.by = NULL, subsample = TRUE, sandwich = TRUE, fix = FALSE, ...)
{
  ## compute meat of sandwich
  rval <- meatPC(x, cluster = cluster, order.by = order.by, subsample = subsample, ...)
    
  ## full sandwich
  if (sandwich) rval <- sandwich(x, meat. = rval)

  ## check (and fix) if sandwich is not positive semi-definite
  if (fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}

meatPC <- function(x, cluster = NULL, order.by = NULL, subsample = TRUE, kronecker = TRUE, ...)
{
  ## extract estimating functions / aka scores
  if (is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  ef <- estfun(x, ...)

  ## FIXME: The return value only has to be set up this way if 
  ## it is filled componentwise. But at the moment the whole matrix is
  ## computed in one big matrix computation.
  ## rval <- matrix(0, nrow = NCOL(ef), ncol = NCOL(ef),
  ##   dimnames = list(colnames(ef), colnames(ef)))

  ## cluster can either be supplied explicitly or
  ## be an attribute of the model  
  if(is.null(cluster)) cluster <- attr(x, "cluster")
  if(is.null(order.by)) order.by <- attr(x, "order.by")

  ## model matrix
  X <- model.matrix(x)
  if (any(alias <- is.na(coef(x)))) X <- X[, !alias, drop = FALSE]
  attr(X, "assign") <- NULL
    
  ## working residuals
  res <- rowMeans(ef/X, na.rm = TRUE)
  res[apply(abs(ef) < .Machine$double.eps, 1L, all)] <- 0

    ## balanced panel
    bal <- split(order.by, cluster)
    balanced <- isTRUE(length(unique(sapply(bal, length))) == 1L)
    if(balanced) {
    ## split residuals by cluster
    e <- split(res, cluster)    
    ## bind into matrix
    e <- do.call("cbind", e)
    ## set up omega
    N <- prod(dim(e))
    t <- dim(e)[1L]
    n <- dim(e)[2L]
    sigma <- crossprod(e) / t
    if(kronecker) {
        omega <- kronecker(sigma, diag(1L, t))           
    } else {
        xx <- split(X, cluster)
        xxCL <- lapply(1L:n, function(i) matrix(xx[[i]], t, ncol(X)))
        omega <- Reduce("+", lapply(1L:n, function(j) Reduce("+", lapply(1L:n, function(i) sigma[j,i] * t(xxCL[[j]]) %*% xxCL[[i]] / N))))
    }
    } 
   if(!balanced) { 
  ## extract balanced subsample
    if(subsample) {
    ## split residuals by cluster
    e <- split(res, cluster)    
    ## cluster sizes
    ne <- sapply(e, length)
    ## extract "full" clusters
    e <- e[ne == max(ne)]
    ## bind into matrix
    e <- do.call("cbind", e)
    ## set up omega
    N <- prod(dim(e))
    t <- dim(e)[1L]
    n <- dim(e)[2L]
    sigma <- crossprod(e) / t
    if(kronecker) {
        X <- X[cluster %in% colnames(e), ]
        omega <- kronecker(sigma, diag(1L, t))           
    } else {
        xx <- split(X, cluster)
        xx <- xx[names(xx) %in% colnames(e)]
        xxCL <- lapply(1L:n, function(i) matrix(xx[[i]], t, ncol(X)))
        omega <- Reduce("+", lapply(1L:n, function(j) Reduce("+", lapply(1L:n, function(i) sigma[j,i] * t(xxCL[[j]]) %*% xxCL[[i]] / N))))
    }
    } else {
    ## use pairwise calculation for omega (Bailey and Katz, 2011)
    N <- dim(X)[1L]    
    pair <- data.frame(cluster = cluster, order.by = order.by)
    pair$res <- res
    pair <- data.frame(pair, X)
    Tij <- expand.grid(cluster = unique(cluster), order.by = unique(order.by))
    pair <- merge(pair, Tij, by = c("cluster", "order.by"), all = TRUE)
    e <- pair[, 3L]
    X <- pair[, 4L:dim(pair)[2L]]
   
    narows <- apply(!is.na.data.frame(X), 1, prod)
    namatrix <- matrix(narows, length(unique(cluster)), length(unique(order.by)), byrow = TRUE)
    namatrix <- t(namatrix)
    denoma <- crossprod(namatrix)
    X[is.na(X)] <- 0L
    
    ## set up omega
    e <- matrix(e, length(unique(cluster)), length(unique(order.by)), byrow = TRUE)
    e[is.na(e)] <- 0L
    e <- t(e)
    sigma <- crossprod(e) / denoma
    if(kronecker) {
        omega <- kronecker(sigma, diag(1L, NROW(e)))
        X <- as.matrix(X)
    } else {
        t <- length(unique(pair$order.by))
        n <- length(unique(pair$cluster))
        xx <- split(X, pair$cluster)
        xxCL <- lapply(1L:n, function(i) as.matrix(xx[[i]], t, ncol(X)))
        omega <- Reduce("+", lapply(1L:n, function(j) Reduce("+", lapply(1L:n, function(i) sigma[j,i] * t(xxCL[[j]]) %*% xxCL[[i]] / N))))
    }
    }
       }
    if(kronecker) {
        rval <- t(X) %*% omega %*% X / N
    } else {
        rval <- omega
    }
  return(rval = rval)
}
