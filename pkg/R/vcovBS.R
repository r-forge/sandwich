.vcovBSenv <- new.env()

vcovBS <- function(x, ...) {
  UseMethod("vcovBS")
}

vcovBS.default <- function(x, cluster = NULL, R = 250, start = FALSE, ..., use = "pairwise.complete.obs")
{
  ## cluster can either be supplied explicitly or
  ## be an attribute of the model...FIXME: other specifications?
  if (is.null(cluster)) cluster <- attr(x, "cluster")

  ## resort to cross-section if no clusters are supplied
  if(is.null(cluster)) {
    cluster <- try(1L:NROW(estfun(x)), silent = TRUE)
    if(inherits(cluster, "try-error")) cluster <- 1L:NROW(model.frame(x))
  }

  ## process formula 'cluster' specification
  if(inherits(cluster, "formula")) {
    cluster_tmp <- expand.model.frame(x, cluster, na.expand = FALSE)
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)[[1L]]
  }
  
  ## handle omitted or excluded observations
  if((n != NROW(cluster)) && !is.null(x$na.action) && (class(x$na.action) %in% c("exclude", "omit"))) {
    cluster <- cluster[-x$na.action]
  }

  ## split clusters
  cl <- split(seq_along(cluster), cluster)
  
  ## set up coefficient matrix
  cf <- coef(x)
  cf <- matrix(rep.int(NA_real_, length(cf) * R), ncol = length(cf),
    dimnames = list(NULL, names(cf)))

  ## use starting values?
  assign(".vcovBSstart", if(isTRUE(start)) coef(x) else NULL, envir = .vcovBSenv)

  ## update on bootstrap samples
  for(i in 1:R) {
      subset <- unlist(cl[sample(names(cl), length(cl), replace = TRUE)])
      assign(".vcovBSsubset", subset, envir = .vcovBSenv)
      up <- if(is.null(.vcovBSenv$.vcovBSstart)) {
        update(x, subset = .vcovBSenv$.vcovBSsubset, ..., evaluate = FALSE)
      } else {
        update(x, subset = .vcovBSenv$.vcovBSsubset, start = .vcovBSenv$.vcovBSstart, ..., evaluate = FALSE)      
      }
      env <- try(environment(terms(x)))
      if(inherits(env, "try-error")) env <- NULL
      up <- eval(up, envir = env, enclos = parent.frame())
      remove(".vcovBSsubset", envir = .vcovBSenv)
      cfi <- coef(up)
      if(is.null(names(cfi))) {
        cf[i, ] <- cfi
      } else {
        cf[i, names(cfi)] <- cfi
      }
  }
  remove(".vcovBSstart", envir = .vcovBSenv)

  # covariance of bootstrap coefficients
  cov(cf, use = use)
}

