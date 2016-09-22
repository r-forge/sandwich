\name{vcovCL}
\alias{vcovCL}
\alias{meatCL}

\title{Clustered Covariance Matrix Estimation}

\description{
  Estimation of one-way and multi-way clustered
  covariance matrices using an object-oriented approach.
}

\usage{
vcovCL(x, cluster = NULL, type = NULL, sandwich = TRUE, fix = FALSE, \dots)
meatCL(x, cluster = NULL, type = NULL, cadjust = TRUE, white = FALSE, \dots)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations
    or a \code{list} (or \code{data.frame}) thereof. By default,
    either \code{attr(x, "cluster")} is used. If that is also \code{NULL}
    each observation is its own cluster.}
  \item{type}{a character string specifying the estimation type (HC0--HC3). 
    The default is to use \code{"HC1"} for \code{"lm"} objects and
    \code{"HC0"} otherwise.}
  \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the meat matrix is returned.}
  \item{fix}{logical. Should the covariance matrix be fixed to be
    positive semi-definite in case it is not?}
  \item{cadjust}{logical. Should a cluster adjustment be applied?}
  \item{white}{logical. Should the HC0 estimate be used for
    the final adjustment in multi-way clustered covariances?}
  \item{\dots}{arguments passed to \code{meatCL}.}
}

\details{
  The function \code{meatCL} is the work horse for estimating
  the meat of clustered sandwich estimators. \code{vcovCL} is a wrapper calling
  \code{\link{sandwich}} and \code{\link{bread}} (Zeileis 2006).
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Petersen MA (2009).
  \dQuote{Estimating Standard Errors in Finance Panel Data Sets: Comparing Approaches},
  \emph{The Review of Financial Studies}, \bold{22}(1), 435--480.
  \doi{10.1093/rfs/hhn053}

Zeileis A (2004).
  \dQuote{Econometric Computing with HC and HAC Covariance Matrix Estimator},
  \emph{Journal of Statistical Software}, \bold{11}(10), 1--17.
  \doi{10.18637/jss.v011.i10}

Zeileis A (2006).
  \dQuote{Object-Oriented Computation of Sandwich Estimators},
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  \doi{10.18637/jss.v016.i09}
}

\seealso{\code{\link{vcovHC}}}

\examples{
## Petersen's data
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

## clustered covariances
## one-way
vcovCL(m, cluster = PetersenCL$firmid)
## one-way with HC2
vcovCL(m, cluster = PetersenCL$firmid, type = "HC2")
## two-way
vcovCL(m, cluster = PetersenCL[, c("firmid", "year")])

## comparison with cross-section sandwiches
## HC0
all.equal(sandwich(m), vcovCL(m, type = "HC0", cadjust = FALSE))
## HC2
all.equal(vcovHC(m, type = "HC2"), vcovCL(m, type = "HC2"))
## HC3
all.equal(vcovHC(m, type = "HC3"), vcovCL(m, type = "HC3"))
}

\keyword{regression}
