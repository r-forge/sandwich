vcovBS <- function(x, cluster = NULL, R = 250, ...)
{
  ## cluster variable (default 1:n)
  if(is.null(cluster)) {
    cluster <- try(1L:NROW(estfun(x)), silent = TRUE)
    if(inherits(cluster, "try-error")) cluster <- 1L:NROW(model.frame(x))
  }
  cl <- split(seq_along(cluster), cluster)
  
  ## set up coefficient matrix
  cf <- coef(x)
  cf <- matrix(rep.int(0, length(cf) * R), ncol = length(cf))

  ## update on bootstrap samples
  for(i in 1:R) {
    ii <- unlist(cl[sample(names(cl), length(cl), replace = TRUE)])
    cf[i, ] <- coef(update(x, subset = ii))
  }
  
  ## covariance of bootstrap coefficients
  cov(cf)
}
