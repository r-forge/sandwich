\name{vcovBS.lm}
\alias{vcovBS.lm}

\title{(Clustered) Bootstrap Covariance Matrix Estimation for lm objects}

\description{
  Estimation of one-way and multi-way clustered
  covariance matrices using bootstrap techniques.
}

\usage{
\method{vcovBS}{default}(x, cluster = NULL, R = 250, multi0 = TRUE, fix = FALSE,
boot_type = "xy", ..., cl = NULL, parallel = "no",
\dots, use = "pairwise.complete.obs")
}

\arguments{
  \item{x}{a fitted model object.}
  \item{cluster}{a variable indicating the clustering of observations
    a \code{list} (or \code{data.frame}) thereof, or a formula specifying which variables from
    the fitted model should be used (see Examples). By default,
    \code{attr(x, "cluster")} is used.}
  \item{R}{integer. Number of bootstrap replications.}
  \item{fix}{logical. Should the covariance matrix be fixed to be
    positive semi-definite in case it is not?}
  \item{boot_type}{a character string specifying the type of bootstrap to use. One of
    "xy", "residual", "wild", "rademacher", "mammen", "norm", "wild-rademacher", 
    "wild-mammen", or "wild-norm". See details for a description of each.}
  \item{multi0}{logical. Should the HC0 estimate be used for
    the final adjustment in multi-way clustered covariances?}
  \item{use}{character. Specification passed to \code{\link[stats]{cov}} for
    handling missing coefficients/parameters.}
  \item{cl}{parallel cluster to be passed to \code{\link{boot}}. See that function for details.}
  \item{parallel}{parallel commands to be passed to \code{\link{boot}}. See that function for details.}
  \item{\dots}{arguments passed to \code{meatCL}.}
}

\details{
  Clustered sandwich estimators are used to adjust inference when errors
  are correlated within (but not between) clusters. See the documentation for \code{\link{vcovCL}} 
  for specifics about standard error clustering. This function allows
  for clustering in arbitrary many cluster dimensions (e.g., firm, time, industry), given all
  dimensions have enough clusters (for more details, see Cameron et al. 2011).
  Unlike \code{vcovCL}, \code{vcovBS.lm} uses a bootstrap rather than an asymptotic solution.
  
  Specifying \code{boot_type = "xy"}, which is the default, implements an "xy" aka "pairs" cluster
  bootstrap. The xy cluster bootstrap resamples with replacement from the observations used in the 
  original model, by cluster instead of resampling by individual observations as in a regular xy 
  bootstrap. For each resampling, new coefficients are estimated via \code{lm.fit()}, the
  collection of which are used to estimate the variance-covariance the function returns.
  
  Specifying \code{boot_type = "residual"} implements a residual cluster bootstrap. The residual
  cluster bootstrap resamples the residuals (as above, by cluster) rather than the observations.
  The resampled residuals are then used to re-form the response variable \code{y* = \hat{y} + e*}
  where \hat{y} is the original fitted values and e* is the resampled residuals. Coefficients are
  reestimated using \code{qr.coef()}. As Cameron et al. (2008) point out, the residual cluster 
  bootstrap is not well-defined when the clusters are unbalanced--different numbers of observations 
  (cluster sizes) among the clusters make it unclear how to correctly apply the residuals of one 
  cluster to another. The function produces a warning if it detects unbalanced clusters.
  
  Specifying \code{boot_type = "wild"} implements a wild cluster bootstrap using the Rademacher
  distribution. A wild bootstrap does not resample anything, but instead reforms the
  dependent variable by multiplying the residual by a randomly drawn value and adding the
  result to the fitted value. There are three built-in distributions to draw multipliers from for wild bootstraps:
  the Rademacher (\code{boot_type = "rademacher"}, the default), which draws from [-1, 1],
  each with p = 0.5, Mammen's suggested distribution (\code{boot_type = "mammen"}, see 
  Mammen, 1993), and the standard normal/Gaussian distribution (\code{boot_type = "norm"}).
  The default is the Rademacher distribution, following Cameron et al. (2008).  Alternatively, you can
  set the function to draw multipliers from by assigning
  \code{boot_type} to a function that takes one argument \code{n}--the number of random values to produce--and
  returns a vector of random values of length \code{n}.
  
  For coding clarity, \code{"wild-rademacher"}, \code{"wild-mammen"}, and \code{"wild-norm"} can be 
  specified as the \code{boot_type} instead; the behavior is identical to specifying (respectively) 
  \code{"rademacher"}, \code{"mammen"}, or \code{"norm"}.

  Cameron et al. (2011) observe that sometimes the covariance matrix is
  not positive-semidefinite. To force the covariance matrix to be
  positive-semidefinite, set
  \code{fix = TRUE}. Following Cameron et al. (2011), the eigendecomposition of the estimated
  covariance matrix is used and any negative eigenvalue(s) are converted to zero.

  
}

\value{
  A matrix containing the covariance matrix estimate.
}

\references{
Cameron AC, Gelbach JB, Miller DL (2008).
  \dQuote{Bootstrap-Based Improvements for Inference with Clustered Errors},
  \emph{The Review of Economics and Statistics}, \bold{90}(3),
  414--427.
  \doi{10.3386/t0344}
  
Cameron AC, Gelbach JB, Miller DL (2011).
  \dQuote{Robust Inference With Multiway Clustering},
  \emph{Journal of Business & Ecomomic Statistics}, \bold{29}(2),
  238--249.
  \doi{10.1198/jbes.2010.07136}
  
Mammen, E. (1993). 
  \dQuote{Bootstrap and wild bootstrap for high dimensional linear models.}, 
  \emph{The Annals of Statistics}, \bold{21}(1),
  255--285.
  \doi{10.1214/aos/1176349025}
}

\seealso{\code{\link{vcovHC}}}

\examples{
\dontrun{
## Petersen's data
data("PetersenCL", package = "sandwich")
fm <- lm(y ~ x, data = PetersenCL)

## clustered covariances
## one-way xy cluster bootstrap
vcovBS(m, cluster = ~ firm)
## two-way wild cluster bootstrap with Mammen distribution
vcovBS(m, cluster = ~ firm + year, boot_type = "wild-mammen")

}
}

\keyword{regression}
\keyword{bootstrap}