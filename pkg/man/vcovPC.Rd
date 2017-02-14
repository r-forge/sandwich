\name{vcovPC}
\alias{vcovPC}
\alias{meatPC}


\title{Panel-Corrected Covariance Matrix Estimation}

\description{
  Estimation of sandwich covariances a la Beck and Katz (1995) for panel data.
}

\usage{
vcovPC(x, cluster = NULL, order.by = NULL,
  subsample = TRUE, sandwich = TRUE, fix = FALSE, \dots)
meatPC(x, cluster = NULL, order.by = NULL,
  subsample = TRUE, kecker = TRUE, \dots)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations
    or a \code{list} (or \code{data.frame}) thereof. By default,
    either \code{attr(x, "cluster")} is used. If that is also \code{NULL}
    each observation is its own cluster.}
  \item{order.by}{a variable indicating the aggregation within time periods.}
  \item{subsample}{logical. FIXME.}
  \item{sandwich}{logical. FIXME.}
  \item{kecker}{logical. FIXME.}
  \item{\dots}{arguments passed to the \code{meatPC} or \code{estfun} function,
    respectively.}
}

\details{
  \code{vcovPC} is a function for estimating Beck and Katz (1995)
  panel-corrected covariance matrix.
 
  The function \code{meatPC} is the work horse for estimating
  the meat of Beck and Katz (1995)
  covariance matrix estimators. \code{vcovPC} is a wrapper calling
  \code{\link{sandwich}} and \code{\link{bread}} (Zeileis 2006).
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Beck N & Katz JN (1995).
\dQuote{What To Do (and Not To Do) with Time-Series-Cross-Section Data
  in Comparative Politics},
\emph{American Political Science Review}, \bold{89}(3), 634--647.
URL \url{http://www.jstor.org/stable/2082979}
  
Zeileis A (2004).
  \dQuote{Econometric Computing with HC and HAC Covariance Matrix Estimator},
  \emph{Journal of Statistical Software}, \bold{11}(10), 1--17.
  \doi{10.18637/jss.v011.i10}

Zeileis A (2006).
  \dQuote{Object-Oriented Computation of Sandwich Estimators},
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  \doi{10.18637/jss.v016.i09}
}

\examples{
## Petersen's data
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

## Beck and Katz (1995) standard errors
## balanced panel
sqrt(diag(vcovPC(m, cluster = PetersenCL$firmid, order.by = PetersenCL$year)))

## Kronecker-product or lapply
system.time(sqrt(diag(vcovPC(m, cluster = PetersenCL$firmid, order.by = PetersenCL$year, kecker = TRUE))))
system.time(sqrt(diag(vcovPC(m, cluster = PetersenCL$firmid, order.by = PetersenCL$year, kecker = FALSE))))

## unbalanced panel
if(require("pcse")) {
data("aglUn", package = "pcse")
agllm <- lm(growth ~ lagg1 + opengdp + openex + openimp + central +
leftc + inter, data = aglUn)
sqrt(diag(vcovPC(agllm, aglUn$country, aglUn$year, subsample = FALSE)))
}
}

\keyword{regression}
