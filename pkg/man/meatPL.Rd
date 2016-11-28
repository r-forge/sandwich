\name{meatPL}
\alias{meatPL}


\title{Clustered Covariance Matrix Estimation for panel data}

\description{
  Estimation of sandwich covariances for panel data.
}

\usage{
vcovPL(x, cluster = NULL, order.by = NULL,
  kernel = "Bartlett", sandwich = TRUE, fix = FALSE, \dots)
meatPL(x, cluster = NULL, order.by = NULL, kernel = "Bartlett", lag =
"max", bw = NULL, adjust = TRUE, \dots)
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations
    or a \code{list} (or \code{data.frame}) thereof. By default,
    either \code{attr(x, "cluster")} is used. If that is also \code{NULL}
    each observation is its own cluster.}
  \item{order.by}{a variable indicating the aggregation within time periods.}
  \item{kernel}{a character specifying the kernel used. All kernels used
  are described in Andrews (1991).}
  \item{sandwich}{logical. Should the sandwich estimator be computed? If set to FALSE only the
meat matrix is returned.}
  \item{lag}{character or logical, indicating the lag length used. By
  default, the maximal lag length \eqn{T - 1} is used.}
  \item{bw}{numeric. The bandwidth of the kernel (corresponds to the
  truncation lag).}
  \item{adjust}{logical. Should a finite sample adjustment be made? This
  amounts to multiplication with \eqn{n/(n - k)} where \eqn{n} is the
  number of observations and \eqn{k} is the number of estimated parameters.}
  \item{\dots}{arguments passed to the \code{estfun} function.}
}

\details{
  \code{vcovPL} is a function for estimating the Newey-West (1978) and Driscoll and Kraay (1998)
  covariance matrix.
  Driscoll and Kraay (1998) apply a Newey-West type correction to the
  sequence of cross-sectional averages of the moment conditions (see
  Hoechle (2007)).
  
  The function \code{meatPL} is the work horse for estimating
  the meat of clustered sandwich estimators. \code{vcovPL} is a wrapper calling
  \code{\link{sandwich}} and \code{\link{bread}} (Zeileis 2006).

  Default lag length is the maximal lag length \eqn{T - 1}. For
  \code{lag = "NW1994"},
  the lag length is taken from the first step of Newey and West's
  (1994) plug-in procedure.
  For \code{lag = "NW1987"}, the lag length is chosen from the heuristic
  \eqn{floor[T^{(1/4)}]}. More details on lag length selection in Hoechle (2007).
  
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Andrews DWK (1991).
\dQuote{Heteroscedasticity and Autocorrelation Consistent Covariance
  Matrix Estimation},
\emph{Econometrica}, 817--858.

Driscoll JC &  Kraay AC (1998).
\dQuote{Consistent Covariance Matrix Estimation with Spatially Dependent
Panel Data},
\emph{The Review of Economics and Statistics}, \bold{80}(4), 549--560.

Hoechle D (2007).
\dQuote{Robust Standard Errors for Panel Regressions with Cross-Sectional Dependence},
\emph{Stata Journal}, \bold{7}(3), 281--312.

Newey WK & West KD (1987).
\dQuote{Hypothesis Testing with Efficient Method of Moments Estimation},
\emph{International Economic Review}, 777-787.

Newey WK & West KD (1994).
\dQuote{Automatic Lag Selection in Covariance Matrix Estimation},
\emph{The Review of Economic Studies}, \bold{61}(4), 631--653.

Petersen MA (2009).
  \dQuote{Estimating Standard Errors in Finance Panel Data Sets: Comparing Approaches},
  \emph{The Review of Financial Studies}, \bold{22}(1), 435--480.
  \doi{10.1093/rfs/hhn053}

Thompson SB (2011).
  \dQuote{Simple Formulas for Standard Errors That Cluster by Both Firm
    and Time},
  \emph{Journal of Financial Economics}, \bold{99}(1), 1--10.
  \doi{10.1016/j.jfineco.2010.08.016}
  
White H (1980).
  \dQuote{A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity},
  \emph{Econometrica}, 817--838.
  \doi{10.2307/1912934}
  
Zeileis A (2004).
  \dQuote{Econometric Computing with HC and HAC Covariance Matrix Estimator},
  \emph{Journal of Statistical Software}, \bold{11}(10), 1--17.
  \doi{10.18637/jss.v011.i10}

Zeileis A (2006).
  \dQuote{Object-Oriented Computation of Sandwich Estimators},
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  \doi{10.18637/jss.v016.i09}
}

\seealso{\code{\link{vcovPL}}}

\examples{
## Petersen's data
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

## Driscoll and Kraay standard errors
## lag length set to maximum lag length T - 1
sqrt(diag(vcovPL(m, cluster = petersen$firmid, adjust =
FALSE)))

## lag length set to rule of thumb proposed by Hoechle (2007)
sqrt(diag(vcovPL(m, cluster = petersen$firmid, lag = "NW1994")))

## lag length set to T^(1/4) 
sqrt(diag(vcovPL(m, cluster = petersen$firmid, lag = "P2009")))
}
}

\keyword{regression}
