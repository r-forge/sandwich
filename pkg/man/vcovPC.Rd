
\name{vcovPC}
\alias{vcovPC}
\alias{meatPC}

\title{Panel-Corrected Covariance Matrix Estimation}

\description{
  Estimation of sandwich covariances a la Beck and Katz (1995) for panel data.
}

\usage{
vcovPC(x, cluster = NULL, order.by = NULL, subsample = TRUE, sandwich = TRUE, fix = FALSE, \dots)
meatPC(x, cluster = NULL, order.by = NULL, subsample = TRUE, kecker = TRUE, \dots)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations
    or a \code{list} (or \code{data.frame}) thereof. By default,
    either \code{attr(x, "cluster")} is used.}
  \item{order.by}{a variable indicating the aggregation within time periods.}
  \item{subsample}{logical. For unbalanced panels. Indicating whether
  the meat should be estimated pair- or casewise.}
  \item{sandwich}{logical. Should the sandwich estimator be computed?
    If set to \code{FALSE} only the meat matrix is returned.}
  \item{fix}{logical. Should the covariance matrix be fixed to be
    positive semi-definite in case it is not?}  
  \item{kecker}{logical. Calculate the meat via the
  Kronecker-product, shortening the computation time for small
  matrices. For large matrices, set \code{kecker = FALSE}.}
  \item{\dots}{arguments passed to the \code{meatPC} or \code{estfun} function,
    respectively.}
}

\details{
  \code{vcovPC} is a function for estimating Beck and Katz (1995)
  panel-corrected covariance matrix.
 
  The function \code{meatPC} is the work horse for estimating
  the meat of Beck and Katz (1995) covariance matrix estimators.
  \code{vcovPC} is a wrapper calling
  \code{\link{sandwich}} and \code{\link{bread}} (Zeileis 2006).

  Following Bailey and Katz (2011), there are two alternatives to
  estimate the meat for unbalanced panels.
  For \code{subsample == TRUE}, a balanced subset of the panel is used,
  whereas for \code{subsample == FALSE}, a pairwise balanced sample is
  employed.
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Bailey D & Katz JN (2011).
\dQuote{Implementing Panel-Corrected Standard Errors in R: The pcse Package},
\emph{Journal of Statistical Software, Code Snippets}, \bold{42}(1), 1--11.
URL \url{http://www.jstatsoft.org/v42/c01/}
  
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

## calculation of the meat via the Kronecker-product 
system.time(sqrt(diag(vcovPC(m, cluster = PetersenCL$firmid, order.by = PetersenCL$year, kecker = TRUE))))
system.time(sqrt(diag(vcovPC(m, cluster = PetersenCL$firmid, order.by = PetersenCL$year, kecker = FALSE))))

## unbalanced panel
if(requireNamespace("pcse", quietly = TRUE)) {
data("aglUn", package = "pcse")
agllm <- lm(growth ~ lagg1 + opengdp + openex + openimp + central +
leftc + inter, data = aglUn)
sqrt(diag(vcovPC(agllm, aglUn$country, aglUn$year, subsample = FALSE)))
}
}

\keyword{regression}