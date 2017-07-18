.vcovBSenv <- new.env()

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
  cf <- matrix(rep.int(0, length(cf) * R), ncol = length(cf),
               dimnames = list(NULL, names(cf)))

  ## update on bootstrap samples
  for(i in 1:R) {
      subset <- unlist(cl[sample(names(cl), length(cl), replace = TRUE)])
      assign(".vcovBSsubset", subset, envir = .vcovBSenv)
      up <- update(x, subset = .vcovBSenv$.vcovBSsubset, evaluate = FALSE)
      env <- try(environment(terms(x)))
      if(inherits(env, "try-error")) env <- NULL
      up <- eval(up, envir = env, enclos = parent.frame())
      remove(".vcovBSsubset", envir = .vcovBSenv)
      cf[i, ] <- coef(up)
  }
  # covariance of bootstrap coefficients
  cov(cf)
}

