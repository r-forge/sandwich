vcovPL <- function(x, cluster = NULL, order.by = NULL,
  kernel = "Bartlett", sandwich = TRUE, fix = FALSE, ...)
{
  ## compute meat of sandwich
  rval <- meatPL(x, cluster = cluster, order.by = order.by, kernel = kernel, ...)

  ## full sandwich
  if(sandwich) rval <- sandwich(x, meat. = rval)

  ## check (and fix) if sandwich is not positive semi-definite
  if(fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}

meatPL <- function(x, cluster = NULL, order.by = NULL,
  kernel = "Bartlett", lag = "NW1987", bw = NULL, adjust = TRUE, aggregate = TRUE, ...) ## adjust/cadjust?
{
  ## extract estimating functions / aka scores
  if (is.list(x) && !is.null(x$na.action)) class(x$na.action) <- "omit"
  ef <- estfun(x, ...)
  k <- NCOL(ef)
  n <- NROW(ef)

  ## cluster can either be supplied explicitly or
  ## be an attribute of the model...FIXME: other specifications?
  if(is.null(cluster)) cluster <- attr(x, "cluster")

  ## cluster (and potentially order.by) may be a formula
  if(inherits(cluster, "formula")) {
    ## merge if both cluster and order.by are formula: ~ cluster + order.by
    if(inherits(order.by, "formula")) {
      cluster_orderby <- ~ cluster + order.by
      cluster_orderby[[2L]][[3L]] <- order.by[[2L]]
      cluster_orderby[[2L]][[2L]] <- cluster[[2L]]
      cluster <- cluster_orderby
      order.by <- NULL
    }
    ## get variable(s) from expanded model frame
    cluster_tmp <- if("Formula" %in% loadedNamespaces()) { ## FIXME to suppress potential warnings due to | in Formula
      suppressWarnings(expand.model.frame(x, cluster, na.expand = FALSE))
    } else {
      expand.model.frame(x, cluster, na.expand = FALSE)
    }
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)

    ## handle omitted or excluded observations
    if(n != NROW(cluster) && !is.null(x$na.action) && class(x$na.action) %in% c("exclude", "omit")) {
      cluster <- cluster[-x$na.action, , drop = FALSE]
    }
  }

  ## cluster can also be a list with both indexes
  if(is.list(cluster)) {
    if(length(cluster) > 1L & is.null(order.by)) order.by <- cluster[[2L]]
    cluster <- cluster[[1L]]
  }

  ## Newey-West type standard errors if neither cluster nor order.by is specified
  if(is.null(cluster) && is.null(order.by)) cluster <- rep(1, n)

  ## resort to cross-section if no clusters are supplied
  if(is.null(cluster)) cluster <- 1L:n

  ## longitudinal time variable
  if(is.null(order.by)) order.by <- attr(x, "order.by")
  if(is.null(order.by)) {
    ix <- make.unique(as.character(cluster))
    ix <- strsplit(ix, ".", fixed = TRUE)
    ix <- sapply(ix, "[", 2L)
    ix[is.na(ix)] <- "0"
    order.by <- as.integer(ix) + 1L
  }

  ## handle omitted or excluded observations in both cluster and/or order.by
  ## (again, in case the formula interface was not used)
  if(any(n != c(NROW(cluster), NROW(order.by))) && !is.null(x$na.action) && class(x$na.action) %in% c("exclude", "omit")) {
    if(n != NROW(cluster))   cluster <-  cluster[-x$na.action]
    if(n != NROW(order.by)) order.by <- order.by[-x$na.action]
  }

  ## catch NAs in cluster/order.by -> need to be addressed in the model object by the user
  if(anyNA(cluster) || anyNA(order.by)) stop("cannot handle NAs in 'cluster' or 'order.by': either refit the model without the NA observations or impute the NAs")

  ## reorder scores and cluster/time variables
  index <- order(order.by)
  ef <- ef[index, , drop = FALSE]
  order.by <- order.by[index]
  cluster <- cluster[index]

  ## aggregate within time periods?
  if(aggregate) {
    nt <- length(unique(order.by))
    if(nt < n) {
      ef <- if(nt == 1L) {
        matrix(colSums(ef), nrow = 1L, dimnames = list(NULL, colnames(ef)))
      } else {
        apply(ef, 2L, tapply, order.by, sum)
      }
    }
    nt <- NROW(ef)
    nc <- length(unique(cluster))
  } else {
    order.by <- as.integer(factor(order.by))
    nt <- order.by[length(order.by)]
    cluster <- as.integer(factor(cluster))
    nc <- cluster[length(cluster)]
    rownames(ef) <- paste(order.by, cluster, sep = "-")
  }  

  ## lag/bandwidth selection
    if(is.character(lag)) {
        if(lag == "P2009") lag <- "max" 
        switch(match.arg(lag, c("max", "NW1987", "NW1994")),
               "max" = {
                   lag <- nt - 1L
               },
               "NW1987" = {
                   lag <- floor(nt^(1/4))
                       },
               "NW1994" = {
                   lag <- floor(4L * (nt / 100L)^(2/9))
                       }
               )}
  if(is.numeric(lag)) lag <- lag
  if(is.null(lag) & is.null(bw)) lag <- nt - 1L
  if(!is.null(lag) & is.null(bw)) bw <- lag + 1

  ## set up kernel weights up to maximal number of lags
  weights <- kweights(0L:(nt - 1L)/bw, kernel = kernel)
  weights <- weights[1L:max(which(abs(weights) > 0))]

  rval <- 0.5 * crossprod(ef) * weights[1L]
    
  if(length(weights) > 1L) {
    for (ii in 2L:length(weights)) {
      rval <- rval + weights[ii] * if(aggregate) {
        ## Driscoll & Kraay (1998)
        crossprod(ef[1L:(nt - ii + 1L), , drop = FALSE], ef[ii:nt, , drop = FALSE])
      } else {
        ## restricting cross-sectional and cross-serial correlation to zero
        ## assure to use only clusters where both lags in the crossprod exist
        ix1 <- paste(rep(1L:(nt - ii + 1L), each = nc), rep(1L:nc, nt - ii + 1L), sep = "-")
        ix2 <- paste(rep(ii:nt, each = nc), rep(1L:nc, nt - ii + 1L), sep = "-")
        ok <- (ix1 %in% rownames(ef)) & (ix2 %in% rownames(ef))
        crossprod(ef[ix1[ok], , drop = FALSE], ef[ix2[ok], , drop = FALSE])
      }
    }
 }
    
  rval <- rval + t(rval)

  if(adjust) rval <- n/(n - k) * rval

  return(rval/n)
}
