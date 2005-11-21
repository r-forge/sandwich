vcovHC <- function(x, 
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4"),
  omega = NULL, sandwich = TRUE, ...)
{
  rval <- meatHC(x, type = type, omega = omega)
  if(sandwich) rval <- sandwich(x, meat = rval, ...)
  return(rval)
}

meatHC <- function(x, 
  type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4"),
  omega = NULL, ...)
{
  ## extract X
  X <- if(is.matrix(x$x)) x$x else model.matrix(x)
  attr(X, "assign") <- NULL
  n <- nrow(X)
  k <- ncol(X)

  ## get residuals and hat values  
  diaghat <- try(hatvalues(x), silent = TRUE)
  res <- if(attr(terms(x), "intercept") > 0) estfun(x)[,1] else rowMeans(estfun(x)/X, na.rm = TRUE)
  
  ## if omega not specified, set up using type
  if(is.null(omega)) {
    type <- match.arg(type)
    if(type == "HC") type <- "HC0"
    switch(type,
      "const" = { omega <- function(residuals, diaghat, df) rep(1, length(residuals)) * sum(residuals^2)/df },
      "HC0"   = { omega <- function(residuals, diaghat, df) residuals^2 },
      "HC1"   = { omega <- function(residuals, diaghat, df) residuals^2 * length(residuals)/df },
      "HC2"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat) },
      "HC3"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat)^2 },
      "HC4"   = { omega <- function(residuals, diaghat, df) residuals^2 / (1 - diaghat)^pmin(4, length(residuals) * diaghat/as.integer(round(sum(diaghat), digits = 0))) })
  }
  
  ## process omega
  if(is.function(omega)) omega <- omega(res, diaghat, n-k)
  rval <- sqrt(omega) * X
  rval <- crossprod(rval)/n

  return(rval)
}
