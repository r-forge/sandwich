bread.default <- function(x, ...) vcov(x) * nobs(x)

## data generating process (dgp)
dgp <- function(nid = 100L, nround = 5L,
  coef = c(0, 0.85, 0.5, 0.7), rho = 0.5, xrho = 0.5,
  dist = "gaussian", type = "copula", link = NULL, ...)
{
  ## match distribution and type of correlation
  dist <- match.arg(dist, c("gaussian", "poisson", "zip", "hurdle", "beta", "logit", "zerotrunc"))
  type <- match.arg(type, c("copula","copula-ar1", "ranef"))
  
  ## sample size
  n <- nid * nround

  ## experimental design variables
  d <- data.frame(
    id = rep(1L:nid, each = nround),
    round = rep(1L:nround, nid)
  )

  ## subject covariates plus random effect
  make_x <- function(corr) {
      rnorm(nid, mean = 0, sd = sqrt(xrho))[d$id] +
      rnorm(n, mean = 0, sd = sqrt(1 - xrho))
  }
  d$x1 <- make_x(corr)
  d$x2 <- rnorm(nid, mean = 0, sd = 1)[d$id]
  d$x3 <- rnorm(n, mean = 0, sd = 1)
  d$ranef <- if(type == "ranef") rnorm(nid, mean = 0, sd = sqrt(rho/(1 - rho)))[d$id] else 0

  ## draw from a normal copula
  if(type == "copula") {
    nc <- copula::normalCopula(rho, dim = nround)
    rcopula <- function(n) {
      rval <- copula::rCopula(n/nround, nc)
      as.vector(t(rval))
    }
  }

  ## draw from a normal copula with AR(1) structure
  if(type == "copula-ar1") {
    nc <- copula::normalCopula(rho, dim = nround, dispstr = "ar1")
    rcopula <- function(n) {
      rval <- copula::rCopula(n/nround, nc)
      as.vector(t(rval))
    }
  }
        
  ## response distribution and link function
  if(is.character(dist)) {
    switch(dist,
      "gaussian" = {
        dist <- if(type == "ranef") {
	  function(n, mu, ...) rnorm(n, mean = mu, ...)
	} else {
	  function(n, mu, ...) qnorm(rcopula(n), mean = mu, ...)
	}
        if(is.null(link)) link <- "identity"
      }, 
      
      "poisson" = {
        dist <- if(type == "ranef") {
	  function(n, mu, ...) rpois(n, lambda = mu)
	} else {
	  function(n, mu, ...) qpois(rcopula(n), lambda = mu)
	}
        if(is.null(link)) link <- "log"
      },
      
      "zip" = {
        dist <- if(type == "ranef") {
	  function(n, mu, ...) countreg::rzipois(n, lambda = mu, pi = 0.3)
	} else {
	  function(n, mu, ...) countreg::qzipois(rcopula(n), lambda = mu, pi = 0.3)
	}
        if(is.null(link)) link <- "log"
      },

      "hurdle" = {
        dist <- if(type == "ranef") {
	  function(n, mu, ...) countreg::rhpois(n, lambda = mu, pi = 0.7)
	} else {
	  function(n, mu, ...) countreg::qhpois(rcopula(n), lambda = mu, pi = 0.7)
	}
        if(is.null(link)) link <- "log"
      },

      "beta" = {
        dist <- if(type == "ranef") {
	  function(n, mu, ...) rbeta(n, shape1 = mu * 10, shape2 = (1 - mu) * 10)
	} else {
	  function(n, mu, ...) qbeta(rcopula(n), shape1 = mu * 10, shape2 = (1 - mu) * 10)
	}
        if(is.null(link)) link <- "logit"
      },

      "logit" = {
        dist <- if(type == "ranef") {
	  function(n, mu, ...) rbinom(n, size = 1, prob = mu)
	} else {
	  function(n, mu, ...) qbinom(rcopula(n), size = 1, prob = mu)
	}
        if(is.null(link)) link <- "logit"
      },

      "zerotrunc" = {
        dist <- if(type == "ranef") {
	  function(n, mu, ...) countreg::rztpois(n, lambda = mu)
	} else {
	  function(n, mu, ...) countreg::qztpois(rcopula(n), lambda = mu)
	}
        if(is.null(link)) link <- "log"
      }

  )}

  if(is.null(link)) link <- "identity"
  if(is.character(link)) link <- make.link(link)
  if(inherits(link, "link-glm")) link <- link$linkinv

  ## compute linear predictor, expectation mu, and response
  d$eta <- coef[1L] + coef[2L] * d$x1 + coef[3L] * d$x2 + coef[4L] * d$x3 + d$ranef
  d$mu <- link(d$eta)
  d$response <- dist(n, d$mu, ...)
  
  ## categorical variables
  d <- transform(d,
    id = factor(id),
    round = factor(round)
  )
  
  ## store coefficients for future reference
  names(coef) <- c("(Intercept)", "x1", "x2", "x3")
  attr(d, "coef") <- coef
  
  return(d)
}
 

## model fitting and covariances
fit <- function(data,
  formula = response ~ x1 + x2 + x3,
  dist = c("gaussian", "poisson", "zip", "hurdle", "beta", "logit", "zerotrunc"),
  vcov = c("without", "HC0", "HC1", "HC2", "HC3", "HC0-id", "HC1-id", "HC2-id", "HC3-id", "fixed", "random", "gee", "dk", "bk"),
  level = 0.95)
{
  
  ## response distributions and vcov types
  dist <- match.arg(dist)
  
  ## pooled model
  m <- switch(dist,
    "gaussian" = lm(formula, data = data),
    "poisson" = glm(formula, data = data, family = poisson),
    "zip" = countreg::zeroinfl(formula, data = data),
    "hurdle" = pscl::hurdle(formula, data = data),
    "beta" = betareg::betareg(formula, data = data, phi = FALSE),
    "logit" = glm(formula, data = data, family = binomial),
    "zerotrunc" = glm(formula, data = data, family = ztpoisson)
    #"zerotrunc" = zerotrunc(formula, data = data, dist = "poisson")
  )
  ## fixed effects
  if("fixed" %in% vcov) {
    formula_fe <- update(formula, . ~ . + id)
    m_fe <- switch(dist,
      "gaussian" = lm(formula_fe, data = data),
      "poisson" = glm(formula_fe, data = data, family = poisson),
      "zip" = countreg::zeroinfl(formula_fe, data = data),
      "hurdle" = pscl::hurdle(formula_fe, data = data),
      "beta" = betareg::betareg(formula_fe, data = data, phi = FALSE),
      "logit" = glm(formula_fe, data, family = binomial),
      "zerotrunc" = glm(formula_fe, data = data, family = ztpoisson)
      #"zerotrunc" = zerotrunc(formula_fe, data = data, dist = "poisson") 
    )
  } else {
    m_fe <- NULL
  }
  ## random effects
  if("random" %in% vcov) {
    formula_re <- update(formula, . ~ . + (1 | id))
    m_re <- switch(dist,
      "gaussian" = lme4::lmer(formula_re, data = data, REML = FALSE),
      "poisson" = lme4::glmer(formula_re, data = data, family = poisson),
      "zip" = NULL,
      "hurdle" = NULL,
      "beta" = NULL,
      "logit" = lme4::glmer(formula_re, data = data, family = binomial),
      "zerotrunc" = NULL
    )
  } else {
    m_re <- NULL
  }
  ## GEE
  if("gee" %in% vcov) {
    m_gee <- switch(dist,
      "gaussian" = geepack::geeglm(formula, data = data, id = id, corstr = "exchangeable", family = gaussian),
      "poisson" = geepack::geeglm(formula, data = data, id = id, corstr = "exchangeable", family = poisson),
      "zip" = NULL,
      "hurdle" = NULL,
      "beta" = NULL,
      "logit" = geepack::geeglm(formula, data = data, id = id, corstr = "exchangeable", family = binomial("logit")),
      "zerotrunc" = NULL
    )
  } else {
    m_gee <- NULL
  }
    
  ## return value: collect coefficients and standard errors
  rval <- data.frame(coef = numeric(0), se = numeric(0), par = character(0),
      vcov = character(0), stringsAsFactors = FALSE)

  if("without" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcov(m))), par = names(coef(m)),
      vcov = "without", stringsAsFactors = FALSE))
  }
  if("HC0" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(sandwich(m))), par = names(coef(m)),
      vcov = "HC0", stringsAsFactors = FALSE))
  }
  if("HC1" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovHC(m, type = "HC1"))), par = names(coef(m)),
      vcov = "HC1", stringsAsFactors = FALSE))
  }
  if("HC2" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovHC(m, type = "HC2"))), par = names(coef(m)),
      vcov = "HC2", stringsAsFactors = FALSE))
  }
  if("HC3" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovHC(m, type = "HC3"))), par = names(coef(m)),
      vcov = "HC3", stringsAsFactors = FALSE))
  }
  if("HC0-id" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovCL(m, cluster = data$id, type = "HC0"))), par = names(coef(m)),
      vcov = "HC0-id", stringsAsFactors = FALSE))
  }
  if("HC1-id" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovCL(m, cluster = data$id, type = "HC1"))), par = names(coef(m)),
      vcov = "HC1-id", stringsAsFactors = FALSE))
  }
  if("HC2-id" %in% vcov & dist != "zip") {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovCL(m, cluster = data$id, type = "HC2"))), par = names(coef(m)),
      vcov = "HC2-id", stringsAsFactors = FALSE))
  }
  if("HC3-id" %in% vcov & dist != "zip") {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovCL(m, cluster = data$id, type = "HC3"))), par = names(coef(m)),
      vcov = "HC3-id", stringsAsFactors = FALSE))
  }
  if("fixed" %in% vcov) {
    k <- length(coef(m))
    rval <- rbind(rval, data.frame(
      coef = coef(m_fe)[1L:k], se = sqrt(diag(vcov(m_fe)))[1L:k], par = names(coef(m_fe))[1L:k],
      vcov = "fixed", stringsAsFactors = FALSE))
  }
  if("random" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = fixef(m_re), se = sqrt(diag(vcov(m_re))), par = names(fixef(m_re)),
      vcov = "random", stringsAsFactors = FALSE))
  }
  if("gee" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m_gee), se = sqrt(diag(m_gee$geese$vbeta)), par = names(coef(m_gee)),
      vcov = "gee", stringsAsFactors = FALSE))
  }

  if("dk" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovPL(m, cluster = data$id, lag = 1, adjust = FALSE))), par = names(coef(m)),
      vcov = "dk", stringsAsFactors = FALSE))
  }

  if("bk" %in% vcov) {
    rval <- rbind(rval, data.frame(
      coef = coef(m), se = sqrt(diag(vcovPC(m, cluster = data$id, order.by = data$round, kecker = TRUE))), par = names(coef(m)),
      vcov = "bk", stringsAsFactors = FALSE))
  }
    
  ## reorder columns
  rownames(rval) <- NULL
  rval <- rval[, c(4, 3, 1, 2)]

  ## FIXME
  rval$par <- gsub("count_", "", rval$par, fixed = TRUE)

  ## further outcomes
  cr <- qnorm((1 - level)/2, lower.tail = FALSE)
  cf <- attr(data, "coef")[rval$par]
  rval$bias <- rval$coef - cf
  rval$mad <- abs(rval$coef - cf)
  rval$power <- as.numeric(abs(rval$coef/rval$se) > cr)
  rval$coverage <- as.numeric(abs(cf - rval$coef)/rval$se < cr)
  return(rval)
}


sim <- function(nrep = 1000, nid = 100L, nround = 5L,
  dist = "gaussian", rho = 0.5, xrho = 0.5,
  coef = c(0, 0.85, 0.5, 0.7), formula = response ~ x1 + x2 + x3,
  vcov = c("without", "HC0", "HC1", "HC2", "HC3", "HC0-id", "HC1-id", "HC2-id", "HC3-id", "fixed", "random", "gee", "dk", "bk"),
  ...,
  cores = NULL)
{
  ## parallelization support
  applyfun <- if(is.null(cores)) {
    lapply
  } else {
    function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = cores)
  }

  ## all factorial combinations of experimental conditions
  par <- expand.grid(nid = nid, nround = nround, dist = dist, rho = rho, xrho = xrho,
    stringsAsFactors = FALSE)

  ## conduct all simulations
  rval <- lapply(1L:nrow(par), function(i) {
    rvali <- applyfun(1L:nrep, function(j) {
      d <- dgp(nid = par$nid[i], nround = par$nround[i], dist = par$dist[i],
        rho = par$rho[i], xrho = par$xrho[i], coef = coef, ...)
      try(fit(d, formula = formula, dist = par$dist[i], vcov = vcov))
    })
    rvali <- rvali[sapply(rvali, class) == "data.frame"]
    rvali[[1L]][, -(1L:2L)] <- Reduce("+", lapply(rvali, "[", , -(1:2)))/length(rvali)
    rvali <- rvali[[1L]]
    rvali$nid <- par$nid[i]
    rvali$nround <- par$nround[i]
    rvali$dist <- par$dist[i]
    rvali$rho <- par$rho[i]
    rvali$xrho <- par$xrho[i]
    return(rvali)
  })
  rval <- do.call("rbind", rval)

  ## turn all experimental condition variables into factors
  rval$dist <- factor(rval$dist)
  rval$vcov <- factor(rval$vcov)
  rval$par <- factor(rval$par)
  rval$nid <- factor(rval$nid)
  rval$nround <- factor(rval$nround)
  rval$rho <- factor(rval$rho)
  rval$xrho <- factor(rval$xrho)

  return(rval)
}

if(FALSE) {

## sandwich package and new vcovCL code
library("sandwich")
library("lme4")
library("geepack")
library("lattice")
library("countreg")
library("VGAM")
source("meatCL.R")
source("simCL.R")
source("vcovPL.R")
source("vcovPC.R")
    
trellis.device(color = FALSE, height = 5, width = 15)

## Experiment 1: How do the different model/covariance estimators
## perform for different types of regressors in the linear regression model?
## Experimental factor: Increase xrho = 0, 0.1, ..., 1
set.seed(1)
s1 <- sim(nrep = 1000, xrho = seq(0, 0.9, by = 0.1), cores = 2, rho = 0.25,
  vcov = c("without", "HC1", "HC1-id", "random", "gee"))

panel.xyref <- function(x, y, ...) {
  panel.abline(h = 0.95, col = 2)
  panel.xyplot(x, y, ...)  
}

xyplot(coverage ~ xrho | par, groups = ~ factor(vcov),
  data = s1, subset = par != "(Intercept)",
  type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)
xyplot(bias ~ xrho | par, groups = ~ factor(vcov),
  data = s1, subset = par != "(Intercept)",
  type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 2: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase xrho = 0, 0.1, ..., 1
set.seed(2)
s2 <- sim(nrep = 1000, coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("gaussian", "poisson", "zip"), vcov = c("without", "HC0", "HC0-id"),
  xrho = seq(0, 1, by = 0.1), cores = 2)

xyplot(coverage ~ xrho | dist, groups = ~ vcov, data = s2,
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)

## Experiment 2a: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase xrho = 0, 0.1, ..., 1
## Increase nid = 400
s2a <- sim(nrep = 1000, nid = 400, coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("gaussian", "poisson", "zip"), vcov = c("without", "HC0", "HC0-id"),
  xrho = seq(0, 1, by = 0.1), cores = 2)

xyplot(coverage ~ xrho | dist + nid, groups = ~ vcov, data = rbind(s2, s2a),
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)

## Experiment 3: How do the different model/covariance estimators
## perform for different types of regressors in the linear regression model?
## Experimental factor: increase nid = 10, 20, ..., 50
set.seed(3)
s3 <- sim(nrep = 1000, nid = seq(10, 50, by = 10), cores = 2,
  vcov = c("without", "HC1", "HC1-id", "fixed", "random", "gee"))

xyplot(coverage ~ nid | par, groups = ~ factor(vcov),
  data = s3, subset = par != "(Intercept)",
  type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 4: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase nid = 10, 20, ..., 50
set.seed(4)
s4 <- sim(nrep = 1000, nid = seq(10, 50, by = 10), coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("gaussian", "poisson", "zip"), vcov = c("without", "HC0", "HC0-id"),
  cores = 2)

xyplot(coverage ~ nid | dist, groups = ~ vcov, data = s4,
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)

## Experiment 5: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase rho = 0, 0.1, ..., 0.9
set.seed(5)
s5 <- sim(nrep = 1000, coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("gaussian", "poisson", "zip"), vcov = c("without", "HC0", "HC0-id"),
  rho = seq(0, 0.9, by = 0.1), cores = 2)

xyplot(coverage ~ rho | dist, groups = ~ vcov, data = s5,
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)

## Experiment 5a: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase rho = 0, 0.1, ..., 0.9
## Increase nid = 400
set.seed(5)
s5a <- sim(nrep = 1000, nid = 400, coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("gaussian", "poisson", "zip"), vcov = c("without", "HC0", "HC0-id"),
  rho = seq(0, 0.9, by = 0.1), cores = 2)

xyplot(coverage ~ rho | dist, groups = ~ vcov, data = s5a,
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)

## Experiment 6: How do the different model/covariance estimators
## perform for different types of regressors in the linear regression model?
## Experimental factor: Increase nround = 5, 10, ..., 25
set.seed(6)
s6 <- sim(nrep = 1000, nround = seq(5, 25, by = 5),
  vcov = c("HC1", "HC2", "HC3", "HC1-id", "HC2-id", "HC3-id"), cores = 2)

xyplot(coverage ~ nround | par, groups = ~ factor(vcov),
  data = s6, subset = par != "(Intercept)",
  type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 7: How do the different model/covariance estimators
## perform for different types of regressors in the linear regression model?
## Experimental factor: increase nround = 5, 10, ..., 25
set.seed(7)
s7 <- sim(nrep = 1000, nround = seq(5, 25, by = 5), cores = 2,
  vcov = c("without", "HC1", "HC1-id", "fixed", "random", "gee"))

xyplot(coverage ~ nround | par, groups = ~ factor(vcov),
  data = s7, subset = par != "(Intercept)",
  type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 8: How do the different model/covariance estimators
## perform for different types of regressors in the linear regression model?
## Experimental factor: Increase rho = 0, 0.1, ..., 0.9
set.seed(8)
s8 <- sim(nrep = 10000,
          vcov = c("without", "HC0", "HC0-id", "random", "gee"),
          rho = seq(0, 0.9, by = 0.1), xrho = 0.25, cores = 2)

xyplot(coverage ~ rho | par, groups = ~ factor(vcov),
  data = s8, subset = par != "(Intercept)",
  type = "b", auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ rho | par, 
  data = s8, subset = par != "(Intercept)",
  type = "b", auto.key = list(columns = 2), panel = panel.xyref)
    

    
## Experiment 9: How do the different covariance estimators
## perform in the Poisson model?
## Experimental factor: increase xrho = 0, 0.1, ..., 0.9
set.seed(9)
s9 <- sim(nrep = 1000, xrho = seq(0, 0.9, by = 0.1), cores = 2, dist = c("poisson"),
          coef = c(0, 0.85, 0, 0), formula = response ~ x1,
          vcov = c("HC1", "HC1-id", "HC2", "HC2-id", "HC3", "HC3-id"))

xyplot(coverage ~ xrho | par, groups = ~ factor(vcov),
       data = s9, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 10: How do the different covariance estimators
## perform in the Poisson model?
## Experimental factor: increase xrho = 0, 0.1, ..., 0.9
set.seed(10)
s10 <- sim(nrep = 1000, rho = seq(0, 0.9, by = 0.1), cores = 2, dist = c("poisson"),
          vcov = c("HC1", "HC1-id", "HC2", "HC2-id", "HC3", "HC3-id"))

xyplot(coverage ~ rho | par, groups = ~ factor(vcov),
       data = s10, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 11: How do the different covariance estimators
## perform in the Poisson model?
## Experimental factor: increase nround = 5, 10, ..., 25
set.seed(11)
s11 <- sim(nrep = 1000, nround = seq(10, 25, by = 5), cores = 2, dist = c("poisson"),
          vcov = c("HC1", "HC1-id", "HC2", "HC2-id", "HC3", "HC3-id"))

xyplot(coverage ~ nround | par, groups = ~ factor(vcov),
       data = s11, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 12: How do the different covariance estimators
## perform in the Gaussian, Poisson, logit and zerotruncation Poisson model?
## Experimental factor: increase nid = 10, 50, 100, ..., 250
set.seed(12)
s12 <- sim(nrep = 10000, nid = c(10, seq(50, 250, by = 50)), cores = 4,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           dist = c("gaussian","poisson", "logit", "zerotrunc"),
           vcov = c("HC0-id","HC1-id","HC2-id","HC3-id"),
           rho = 0.25, xrho = 0.25)

xyplot(coverage ~ nid | dist, groups = ~ factor(vcov),
       data = s12, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ nid | paste(dist, par), groups = ~ factor(vcov),
       data = s12, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)

set.seed(12)
s12a <- sim(nrep = 1000, nid = c(10, seq(50, 250, by = 50)), cores = 4,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           dist = c("gaussian", "poisson", "logit", "zerotrunc"),
           vcov = c("HC0-id","HC1-id","HC2-id","HC3-id"),
           rho = 0.1, xrho = 0.25)
xyplot(coverage ~ nid | dist, groups = ~ factor(vcov),
       data = s12a, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)
    
## Experiment 13: How do the different covariance estimators
## perform in the Poisson model?
## Experimental factor: increase xrho = 0, 0.1, ..., 0.9
set.seed(13)
s13 <- sim(nrep = 1000, xrho = seq(0, 0.9, by = 0.1), cores = 2, dist = c("poisson"),
          vcov = c("HC1", "HC1-id", "HC2", "HC2-id", "HC3", "HC3-id"))

xyplot(coverage ~ xrho | par, groups = ~ factor(vcov),
       data = s13, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 14: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase nid = 50, 60, ..., 100
## nround = 10
set.seed(14)
s14 <- sim(nrep = 1000, nid = seq(50, 100, by = 10), nround = 10, coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("gaussian", "poisson", "zip"), vcov = c("without", "HC0", "HC0-id"),
  cores = 2)

xyplot(coverage ~ nid | dist, groups = ~ vcov, data = s14,
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)

## Experiment 15: How do the different covariance estimators
## perform in the linear model?
## Experimental factor: increase nid = 10, 50, 100, ..., 250
set.seed(15)
s15 <- sim(nrep = 1000, nid = c(10, seq(50, 250, by = 50)), cores = 2, dist = c("gaussian"),
          vcov = c("HC0-id","HC1-id","HC2-id","HC3-id"), rho = 0.25, xrho = 0.25)

xyplot(coverage ~ nid | par, groups = ~ factor(vcov),
       data = s15, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 16: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase rho = 0, 0.1, ..., 0.9
set.seed(16)
s16 <- sim(nrep = 1000, coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("poisson", "zip", "hurdle"), vcov = c("without", "HC0", "HC0-id"),
  rho = seq(0, 0.9, by = 0.1))

xyplot(coverage ~ rho | dist, groups = ~ vcov, data = s16,
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)

## Experiment 17: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase rho = 0, 0.1, ..., 0.9
set.seed(17)
s17 <- sim(nrep = 1000, nid = 400, coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("poisson", "zip", "hurdle"), vcov = c("without", "HC0", "HC0-id"),
  rho = seq(0, 0.9, by = 0.1))

xyplot(coverage ~ rho | dist, groups = ~ vcov, data = s17,
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)


## Experiment 18: How do clustered covariances perform for different
## types of responses?
## Experimental factor: Increase rho = 0, 0.1, ..., 0.9
set.seed(18)
s18 <- sim(nrep = 1000, nid = 800, coef = c(0, 0.85, 0, 0), formula = response ~ x1,
  dist = c("poisson", "zip", "hurdle"), vcov = c("without", "HC0", "HC0-id"),
  rho = seq(0, 0.9, by = 0.1))

xyplot(coverage ~ rho | dist, groups = ~ vcov, data = s18,
  subset = par == "x1", type = "b", auto.key = TRUE, panel = panel.xyref)

## Experiment 19: How do the different model/covariance estimators
## perform for different types of regressors in the beta regression model?
## Experimental factor: Increase rho = 0, 0.1, ..., 1
set.seed(19)
s19 <- sim(nrep = 1000, dist = "beta", vcov = c("without", "HC0", "HC0-id"),
           rho = seq(0, 0.9, by = 0.1), cores = 3)

xyplot(coverage ~ rho | par, groups = ~ factor(vcov),
  data = s19, subset = par != "(Intercept)", type = "b",
  xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 20: How do the different model/covariance estimators
## perform for different types of regressors in the beta regression model?
## Experimental factor: Increase rho = 0, 0.1, ..., 1
set.seed(20)
s20 <- sim(nrep = 1000, dist = "beta", vcov = c("without", "HC0", "HC0-id"),
           rho = seq(0, 0.9, by = 0.1), xrho = 0.25, cores = 3)

xyplot(coverage ~ rho | par, groups = ~ factor(vcov),
  data = s20, subset = par != "(Intercept)", type = "b",
  xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ rho | par,
  data = s20, subset = par != "(Intercept)", type = "b",
  xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)
    
## Experiment 21: How do the different model/covariance estimators
## perform for different types of regressors in the beta regression model?
## Experimental factor: Increase increase nid = 10, 50, 100, ..., 250
set.seed(21)
s21 <- sim(nrep = 1000, dist = "beta", vcov = c("HC0-id", "HC1-id", "HC2-id", "HC3-id"),
           nid = c(10, seq(50, 250, by = 50)))

xyplot(coverage ~ nid | par, groups = ~ factor(vcov),
  data = s21, subset = par != "(Intercept)", type = "b",
  xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 22: How do the different model/covariance estimators
## perform for different types of regressors in the logit model?
## Experimental factor: Increase nid = 10, 50, 100, ..., 250
set.seed(22)
s22 <- sim(nrep = 1000, nid = c(10, seq(50, 250, by = 50)), cores = 2, dist = "logit",
           vcov = c("without", "HC1", "HC1-id", "fixed", "random", "gee"))

xyplot(coverage ~ nid | par, groups = ~ factor(vcov),
  data = s22, subset = par != "(Intercept)",
  type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)
    
## Experiment 23: How do the different model/covariance estimators
## perform for different types of regressors in the logit model?
## Experimental factor: Increase rho = 0, 0.1, ..., 1
set.seed(23)
s23 <- sim(nrep = 1000, dist = "logit", vcov = c("without", "HC0", "HC1-id", "HC3-id"),
           rho = seq(0, 0.9, by = 0.1), cores = 3)

xyplot(coverage ~ rho | par, groups = ~ factor(vcov),
  data = s23, subset = par != "(Intercept)", type = "b",
  xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 24: How do the different model/covariance estimators
## perform for different types of regressors in the logit model?
## Experimental factor: Increase xrho = 0, 0.1, ..., 0.9
set.seed(24)
s24 <- sim(nrep = 1000, dist = "logit", vcov = c("without", "HC0", "HC1-id", "HC3-id"),
           xrho = seq(0, 0.9, by = 0.1), cores = 3)

xyplot(coverage ~ xrho | par, groups = ~ factor(vcov),
  data = s24, subset = par != "(Intercept)", type = "b",
  xlab = expression(paste(rho ["x"])), auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 25: How do the different model/covariance estimators
## perform for different types of regressors in the zero truncated Poisson model?
## Experimental factor: Increase xrho = 0, 0.1, ..., 0.9
set.seed(25)
s25 <- sim(nrep = 1000, dist = "zerotrunc", vcov = c("without", "HC0", "HC1-id", "HC3-id"),
           xrho = seq(0, 0.9, by = 0.1), cores = 3)

xyplot(coverage ~ xrho | par, groups = ~ factor(vcov),
  data = s25, subset = par != "(Intercept)", type = "b",
  xlab = expression(paste(rho ["x"])), auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 26: How do the different model/covariance estimators
## perform for different types of regressors in the beta regression model?
## Experimental factor: Increase xrho = 0, 0.1, ..., 1
set.seed(26)
s26 <- sim(nrep = 1000, dist = "beta", vcov = c("without", "HC0", "HC0-id"),
           xrho = seq(0, 0.9, by = 0.1), cores = 3)

xyplot(coverage ~ xrho | par, groups = ~ factor(vcov),
  data = s26, subset = par != "(Intercept)", type = "b",
  xlab = expression(paste(rho ["x"])), auto.key = list(columns = 2), panel = panel.xyref)


## Experiment 27: How do the different model/covariance estimators
## perform for different types of regressors in the zero truncated Poisson model?
## Experimental factor: Increase rho = 0, 0.1, ..., 0.9
set.seed(27)
s27 <- sim(nrep = 1000, dist = "zerotrunc", vcov = c("without", "HC0", "HC1-id", "HC2-id", "HC3-id"),
           rho = seq(0, 0.9, by = 0.1), xrho = 0.1, cores = 3)

xyplot(coverage ~ rho | par, groups = ~ factor(vcov),
  data = s27, subset = par != "(Intercept)", type = "b",
  xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ rho | par,
  data = s27, subset = par != "(Intercept)", type = "b",
  xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

    
## Experiment 28: How do the different model/covariance estimators
## perform for different types of regressors in the zero truncated Poisson model?
## Experimental factor: Increase nid = 10, 50, 100, ..., 250
set.seed(28)
s28 <- sim(nrep = 1000, dist = "zerotrunc", vcov = c("without", "HC0", "HC1-id", "HC3-id"),
           nid = c(10, seq(50, 250, by = 50)), cores = 3)

xyplot(coverage ~ nid | par, groups = ~ factor(vcov),
  data = s28, subset = par != "(Intercept)", type = "b",
  xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)


## Experiment 29: How do the different model/covariance estimators
## perform for different types of models?
## Experimental factor: Increase rho = 0, 0.1, ..., 0.9. Set xrho to 0.25.
set.seed(29)
s29 <- sim(nrep = 10000, rho = seq(0, 0.9, by = 0.1), xrho = 0.25,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           cores = 3, dist = c("poisson", "logit", "zerotrunc"),
           vcov = c("without", "HC0", "HC0-id"))

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
  data = s29, subset = par != "(Intercept)",
  type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ rho | dist,
  data = s29, subset = par != "(Intercept)",
  type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

    
set.seed(30)
s30 <- sim(nrep = 1000, nid = c(10, seq(50, 250, by = 50)), xrho = 0.1,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           cores = 3, dist = c("poisson", "logit", "zerotrunc"),
           vcov = c("without", "HC1", "HC1-id"))

xyplot(coverage ~ nid | dist, groups = ~ factor(vcov),
  data = s30, subset = par != "(Intercept)",
  type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ nid | dist,
  data = s30, subset = par != "(Intercept)",
  type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)


## Experiment 31: How do the different model/covariance estimators
## perform for different types of models ?
## Experimental factor: Increase rho = 0, 0.1, ..., 1
set.seed(31)
s31 <- sim(nrep = 10000, dist = c("beta", "zip"), vcov = c("without", "HC0", "HC0-id"),
           rho = seq(0, 0.9, by = 0.1), xrho = 0.25, nid = 100, cores = 3,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1)

g1 <- xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
  data = s31, subset = par != "(Intercept)", type = "b", ylab = "coverage",
  xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)
g2 <- xyplot(bias ~ rho | dist,
  data = s31, subset = par != "(Intercept)", type = "b",
  xlab = expression(rho), ylab = "bias", auto.key = list(columns = 2), panel = panel.xyref)    

print(g1, position = c(0, 0, 1, 1),split=c(1,1,1,2), more = TRUE)
print(g2, position = c(0, 0, 1, 1), split = c(1,2,1,2))    

install.packages("latticeExtra")
library("latticeExtra")
doubleYScale(g1, g2, style1 = 0, style2 = 3, add.ylab2 = TRUE, columns = 2, rows = 2)
xyplot(bias + coverage ~ rho | dist, groups = ~ factor(vcov),
  data = s31, subset = par != "(Intercept)", type = "b", ylab = "coverage",
  xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 121: How do the different covariance estimators
## perform in the Gaussian and logit model?
## Experimental factor: increase rho = 0, 0.1, ..., 0.9
set.seed(121)
s121 <- sim(nrep = 10000, nid = 100, nround = 5, cores = 4,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           dist = c("gaussian", "logit2"),
           vcov = c("without", "HC0", "HC0-id"),
           rho = seq(0, 0.9, by = 0.1), xrho = 0.25)

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
data = s121, subset = par != "(Intercept)",
type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ rho | dist, groups = ~ factor(vcov),
data = s121, subset = par != "(Intercept)",
type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment 122: How do the different covariance estimators
## perform in the Gaussian and logit model?
## Experimental factor: increase nid = 10, 50, 100, ..., 250
set.seed(122)
s122 <- sim(nrep = 10000, nid = c(10, seq(50, 250, by = 50)), nround = 5, cores = 4,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           dist = c("gaussian", "logit2"),
           vcov = c("without", "HC0", "HC0-id"),
           rho = 0.25, xrho = 0.25)

xyplot(coverage ~ nid | dist, groups = ~ factor(vcov),
data = s122, subset = par != "(Intercept)",
type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ nid | dist, groups = ~ factor(vcov),
data = s122, subset = par != "(Intercept)",
type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)

    
## Experiment 123: How do the different covariance estimators
## perform in the Gaussian, logit model, and Poisson model?
## Experimental factor: increase rho = 0, 0.1, ..., 0.9
set.seed(123)
s123 <- sim(nrep = 10000, nid = 100, nround = 5,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           dist = c("gaussian", "logit", "poisson"),
           type = c("copula"),
           vcov = c("without", "random", "gee", "HC0", "HC0-id"),
           rho = seq(0, 0.9, by = 0.1), xrho = 0.5)

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
data = s123, subset = par != "(Intercept)",
type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ rho | dist, groups = ~ factor(vcov),
data = s123, subset = par != "(Intercept)",
type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)
    
## Experiment 124: How do the different covariance estimators
## perform in the Gaussian, logit model, and Poisson model?
## Experimental factor: increase rho = 0, 0.1, ..., 0.9
set.seed(124)
s124 <- sim(nrep = 1000, nid = 100, nround = 5,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           dist = c("gaussian", "logit", "poisson"),
           type = c("ranef"),
           vcov = c("without", "gee", "HC0", "HC0-id"),
           rho = seq(0, 0.9, by = 0.1), xrho = 0.5)

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
data = s124, subset = par != "(Intercept)",
type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ rho | dist, groups = ~ factor(vcov),
data = s124, subset = par != "(Intercept)",
type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)

set.seed(125)
s125 <- sim(nrep = 10000, nid = 100, nround = 5,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           dist = c("zip", "zerotrunc"),
           type = c("copula"),
           vcov = c("without", "HC0", "HC0-id"),
           rho = seq(0, 0.9, by = 0.1), xrho = 0.5)    

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
data = s125, subset = par != "(Intercept)",
type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ rho | dist, groups = ~ factor(vcov),
data = s125, subset = par != "(Intercept)",
type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)

set.seed(126)
s126 <- sim(nrep = 10000, nid = c(10, seq(50, 250, by = 50)), cores = 4,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           dist = c("gaussian","poisson", "logit", "zerotrunc"),
           type = "copula",
           vcov = c("HC0-id","HC1-id","HC2-id","HC3-id"),
           rho = 0.25, xrho = 0.25)

xyplot(coverage ~ nid | dist, groups = ~ factor(vcov),
       data = s126, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)

xyplot(bias ~ nid | paste(dist, par), groups = ~ factor(vcov),
       data = s126, subset = par != "(Intercept)",
       type = "b", auto.key = list(columns = 2), panel = panel.xyref)



    
## Experiment Driscoll & Kraay standard errors with vcovPL()
set.seed(127)
s127 <- sim(nrep = 10000, nid = 100, nround = 5, cores = 3,
          coef = c(0, 0.85, 0 , 0), formula = response ~ x1,
          dist = c("gaussian", "poisson", "logit"),
          type = "copula-ar1",
          vcov = c("without", "HC0", "HC0-id", "dk"),
          rho = seq(0, 0.9, by = 0.1), xrho = 0.25)

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
  data = s127, subset = par != "(Intercept)",
  type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)

set.seed(128)
s128 <- sim(nrep = 10000, nid = c(10, seq(50, 250, by = 50)), nround = 5, cores = 3, 
          coef = c(0, 0.85, 0 , 0), formula = response ~ x1,
          dist = c("gaussian", "poisson", "logit"),
          type = "copula-ar1",
          vcov = c("without", "HC0", "HC0-id", "dk"),
          rho = 0.25, xrho = 0.25)

xyplot(coverage ~ nid | dist, groups = ~ factor(vcov),
  data = s128, subset = par != "(Intercept)",
  type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)


set.seed(129)
s129 <- sim(nrep = 1000, rho = seq(0, 0.9, by = 0.1), xrho = 0.25,
           coef = c(0, 0.85, 0, 0), formula = response ~ x1,
           cores = 3, dist = c("poisson", "logit", "zerotrunc"),
           type = "copula",
           vcov = c("without", "HC0", "HC0-id", "dk"))    

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
  data = s129, subset = par != "(Intercept)",
  type = "b", xlab = "G", auto.key = list(columns = 2), panel = panel.xyref)

## Experiment Beck & Katz standard errors with vcovPC()
set.seed(130)
s130 <- sim(nrep = 1000, nid = 100, nround = 5, cores = 3,
          coef = c(0, 0.85, 0 , 0), formula = response ~ x1,
          dist = c("gaussian", "poisson", "logit"),
          type = "copula-ar1",
          vcov = c("without", "HC0", "HC0-id", "dk", "bk"),
          rho = seq(0, 0.9, by = 0.1), xrho = 0.25)

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
  data = s130, subset = par != "(Intercept)",
  type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)


set.seed(131)
s131 <- sim(nrep = 1000, nid = 100, nround = 5, cores = 3,
          coef = c(0, 0.85, 0 , 0), formula = response ~ x1,
          dist = c("gaussian", "poisson", "logit"),
          type = "copula",
          vcov = c("without", "HC0", "HC0-id", "dk", "bk"),
          rho = seq(0, 0.9, by = 0.1), xrho = 0.25)

xyplot(coverage ~ rho | dist, groups = ~ factor(vcov),
  data = s131, subset = par != "(Intercept)",
  type = "b", xlab = expression(rho), auto.key = list(columns = 2), panel = panel.xyref)
    
}
