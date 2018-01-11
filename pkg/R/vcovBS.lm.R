vcovBS.lm.default <- function(x, cluster = NULL, parallel = FALSE, multi0 = NULL, 
                         fix = FALSE, R = 250, boot_type = "xy", wild_type = "rademacher",
                         debug = FALSE, ..., use = "pairwise.complete.obs") {
  
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
  
  n <- NROW(model.frame(x))
  
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
  
  if(is.null(multi0)) {
    if(cluster_dims > 1 && df[tcc, "M"] == prod(df[-tcc, "M"])) {
      multi0 <- TRUE
    } else {
      multi0 <- FALSE
    }
  }
  
  if(multi0) {
    df <- df[-tcc,]
    tcc <- tcc - 1
  }
  
  if(debug) {
    print(acc)
    print(paste("Original Cluster Dimensions", cluster_dims))
    print(paste("Theoretical Cluster Combinations", tcc))
    print(paste("Use White VCOV for final matrix?", multi0))
  }
  
  
  if("model" %in% names(x)) {
    full_data <- x$model
  } else {
    full_data <- model.frame(x)
  }
  
  boot.outs <- list()
  
  par_type <- "no"
  if(length(parallel) > 1) {
    par_cluster <- parallel
    clusterExport(par_cluster, varlist = c("cluster", "x"), envir = environment())
    par_type <- "snow"
  } else if(parallel == TRUE || parallel == "snow") {
    par_type <- "snow"
    par_cluster <- NULL
  } else if(parallel == "multicore") {
    par_type <- "multicore"
    par_cluster <- NULL
  }
  
  dat_loc <- which(names(x$call) == "data") 
  args <- x$call[c(-1, -dat_loc)]
  wild_func <- function() NULL
  
  boot_type <- match.arg(boot_type, c("xy", "residual", "wild"))
  
  # Setup the bootstrap functions, depending upon the type of bootstrap requested
  if(boot_type == "xy") {
    est.func <- cmpfun(function(grp, i, data2, clustvar, reg_arglist, boot_args) {
      j <- unlist(lapply(i, function(n) which(n == clustvar)))
      coef(boot_args$estimator(reg_arglist, data = data2[j,]))
    })
  } else if(boot_type == "residual") {
    args$formula <- update.formula(formula(x), y_boot ~ .)
    full_data$y_boot <- 0
    
    est.func <- cmpfun(function(grp, i, data2, clustvar, reg_arglist, boot_args) {
      j <- unlist(lapply(i, function(n) which(n == clustvar)))
      data2$y_boot <- fitted(boot_args$x) + residuals(boot_args$x)[j]
      coef(boot_args$estimator(reg_arglist, data = data2))
    })
  } else if(boot_type == "wild") {
    args$formula <- update.formula(formula(x), y_boot ~ .)
    full_data$y_boot <- 0
    
    if(is.function(wild_type)) {
      wild_func <- wild_type
    } else {
      wild_type <- match.arg(wild_type, c("rademacher", "mammen", "norm"))
      if(wild_type == "rademacher") {
        wild_func <- cmpfun(function() sample(c(-1, 1), 1))
      } else if(wild_type == "mammen") {
        wild_func <- cmpfun(function() sample(c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2), 1,
                                              prob = c((sqrt(5) + 1) / (2 * sqrt(5)), 
                                                       (sqrt(5) - 1) / (2 * sqrt(5)))))
      } else if(wild_type == "norm") {
        wild_func <- cmpfun(function() rnorm(1))
      }
    }
    
    
    est.func <- cmpfun(function(grp, i, data2, clustvar, reg_arglist, boot_args) {
      j <- unlist(lapply(grp, function(n) rep_len(boot_args$wild_func(), sum(n == clustvar))))
      data2$y_boot <- fitted(x) + residuals(x) * j
      coef(boot_args$estimator(reg_arglist, data = data2))
    })
  }
  boot_args <- new.env()
  boot_args$estimator <- eval(x$call[[1]])
  boot_args$x <- x
  boot_args$wild_func <- wild_func
  
  for(i in 1:tcc) {
    boot.outs[[i]] <- boot(unique(cluster[,i]), est.func, R = R,
                           parallel = par_type, cl = par_cluster,
                           data2 = full_data, clustvar = cluster[,i],
                           reg_arglist = args, boot_args = boot_args)
  }
  
  if(debug) {
    print(df)
    print(vcov_sign)
  }
  
  vcov_matrices <- list()
  for(i in 1:tcc) {
    vcov_matrices[[i]] <- vcov_sign[i] * cov(boot.outs[[i]]$t, use = use)
  }
  
  if(multi0) {
    i <- i + 1
    vcov_matrices[[i]] <- vcov_sign[i] * vcovHC(x, type = "HC0")
  }
  
  if(debug) {
    print(vcov_matrices)
  }
  
  vcov_matrix <- Reduce('+', vcov_matrices)
  
  if(fix) {
    decomp <- eigen(vcov_matrix, symmetric = TRUE)
    if(debug) print(decomp$values)
    pos_eigens <- pmax(decomp$values, rep.int(0, length(decomp$values)))
    vcov_matrix <- decomp$vectors %*% diag(pos_eigens) %*% t(decomp$vectors)
  }
  
  return(vcov_matrix)
}