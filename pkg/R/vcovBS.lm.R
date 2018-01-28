vcovBS.lm <- function(x, cluster = NULL, R = 250, multi0 = TRUE, 
                      fix = FALSE, boot_type = "xy", wild_type = "rademacher",
                      ..., use = "pairwise.complete.obs", 
                      cl = NULL, parallel = "no") {
  
  debug <- FALSE
  arglist <- list(...)
  if("debug" %in% names(arglist)) {
    debug <- arglist$debug
  }
  
  if(inherits(cluster, "formula")) {
    cluster_tmp <- expand.model.frame(x, cluster, na.expand = FALSE)
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
  } else {
    cluster <- as.data.frame(cluster, stringsAsFactors = FALSE)
  }
  
  cluster_dims <- ncol(cluster)
  
  # total cluster combinations, 2^D - 1
  tcc <- 2 ** cluster_dims - 1
  # all cluster combinations
  acc <- list()
  for(i in 1:cluster_dims) {
    acc <- append(acc, combn(1:cluster_dims, i, simplify = FALSE))
  }
  if(debug) print(acc)
  
  # We need to subtract matrices with an even number of combinations and add
  # matrices with an odd number of combinations
  vcov_sign <- sapply(acc, function(i) (-1) ** (length(i) + 1))
  
  # Drop the original cluster vars from the combinations list
  acc <- acc[-1:-cluster_dims]
  if(debug) print(acc)
  
  n <- nobs(x)
  
  # Handle omitted or excluded observations
  if(n != NROW(cluster) && !is.null(x$na.action)) {
    if(class(x$na.action) %in% c("exclude", "omit")) {
      cluster <- cluster[-x$na.action,]
    }
    cluster <- as.data.frame(cluster)  # FIXME: silly error somewhere
  }
  if(debug) print(class(cluster))
  
  # Factors in our clustering variables can potentially cause problems
  # Blunt fix is to force conversion to characters
  i <- !sapply(cluster, is.numeric)
  cluster[i] <- lapply(cluster[i], as.character)
  
  # Make all combinations of cluster dimensions
  if(cluster_dims > 1) {
    for(i in acc) {
      cluster <- cbind(cluster, Reduce(paste0, cluster[,i]))
    }
  }
  
  df <- data.frame(M = integer(tcc),
                   N = integer(tcc),
                   K = integer(tcc))
  
  for(i in 1:tcc) {
    df[i, "M"] <- length(unique(cluster[,i]))
    df[i, "N"] <- length(cluster[,i])
    df[i, "K"] <- x$rank
  }
  
  if(multi0 && tcc > 1) {
    df <- df[-tcc,]
    tcc <- tcc - 1
  }
  
  if(debug) {
    print(acc)
    print(paste("Original Cluster Dimensions", cluster_dims))
    print(paste("Theoretical Cluster Combinations", tcc))
    print(paste("Use White VCOV for final matrix?", multi0))
  }
  
  
  boot.outs <- list()

  boot_type <- match.arg(boot_type, c("xy", "residual", "wild"))
  
  boot_args <- new.env()
  boot_args$x <- x
  boot_args$wild_func <- function() NULL
  
  # Setup the bootstrap functions, depending upon the type of bootstrap requested
  if(boot_type == "xy") {
    boot_args$mx <- model.matrix(x)
    boot_args$y <- model.response(boot_args$x$model)
    est.func <- cmpfun(function(grp, i, clustvar, boot_args) {
      j <- unlist(lapply(i, function(n) which(n == clustvar)))
      coef(lm.fit(x = boot_args$mx[j, ], y = boot_args$y[j]))
    })
  } else if(boot_type == "residual") {
    est.func <- cmpfun(function(grp, i, clustvar, boot_args) {
      j <- unlist(lapply(i, function(n) which(n == clustvar)))
      y_boot <- fitted(boot_args$x) + residuals(boot_args$x)[j]
      qr.coef(boot_args$x$qr, y_boot)
    })
  } else if(boot_type == "wild") {
    
    if(is.function(wild_type)) {
      boot_args$wild_func <- wild_type
    } else {
      wild_type <- match.arg(wild_type, c("rademacher", "mammen", "norm"))
      if(wild_type == "rademacher") {
        boot_args$wild_func <- cmpfun(function(n) sample(c(-1, 1), n, replace = TRUE))
      } else if(wild_type == "mammen") {
        boot_args$wild_func <- cmpfun(function(n) sample(c(-(sqrt(5) - 1) / 2, 
                                                           (sqrt(5) + 1) / 2), 
                                                         n,
                                                         replace = TRUE,
                                                         prob = c((sqrt(5) + 1) / (2 * sqrt(5)), 
                                                                  (sqrt(5) - 1) / (2 * sqrt(5)))))
      } else if(wild_type == "norm") {
        boot_args$wild_func <- cmpfun(function(n) rnorm(n))
      }
    }
    
    clust <- split(seq_along(PetersenCL$year), PetersenCL$year)
    boot_args$grp_indices <- unlist(clust)
    boot_args$grp_sizes <- lapply(clust, length)
    
    est.func <- cmpfun(function(grp, i, clustvar, boot_args) {
      draws <- boot_args$wild_func(length(grp))
      expanded_draws <- unlist(mapply(rep, draws, boot_args$grp_sizes, SIMPLIFY = FALSE))
      y_boot <- fitted(boot_args$x) + residuals(boot_args$x) * expanded_draws[boot_args$grp_indices]
      qr.coef(boot_args$x$qr, y_boot)
    })
  }
  
  for(i in 1:tcc) {
    boot.outs[[i]] <- boot(unique(cluster[,i]), est.func, R = R,
                           parallel = parallel, cl = cl,
                           clustvar = cluster[,i],
                           boot_args = boot_args)
  }
  
  if(debug) {
    print(df)
    print(vcov_sign)
  }
  
  vcov_matrices <- list()
  for(i in 1:tcc) {
    vcov_matrices[[i]] <- vcov_sign[i] * cov(boot.outs[[i]]$t, use = use)
  }
  
  if(multi0 && tcc > 1) {
    i <- i + 1
    vcov_matrices[[i]] <- vcov_sign[i] * sandwich(fm)
  }
  
  if(debug) {
    print(vcov_matrices)
  }
  
  vcov_matrix <- Reduce('+', vcov_matrices)
  
  if(fix && any((eig <- eigen(vcov_matrix, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    if(debug) print(eig$values)
    vcov_matrix <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  
  return(vcov_matrix)
}