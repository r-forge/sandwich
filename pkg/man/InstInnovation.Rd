\name{InstInnovation}
\alias{InstInnovation}

\title{Innovation and Institutional Ownership}

\description{
  Firm-level data on innovation and institutional ownership from 1991 to
  1999 over 803 firms.
}

\usage{data("InstInnovation")}

\format{
A data frame containing 6208 observations on 25 variables.
\describe{
  \item{company}{factor. company names.}
  \item{sales}{numeric. sales (in millions of dollars).}
  \item{acompetition}{numeric. constant inverse Lerner index.}
  \item{competition}{numeric. varying inverse Lerner index.}
  \item{capital}{numeric. net stock of property, plant, and equipment.}
  \item{cites}{integer. future cite-weighted patents.}
  \item{precites}{numeric. presample average of cite-weighted patents.}
  \item{dprecites}{factor. indicates zero precites.}
  \item{patents}{integer. granted patents.}
  \item{drandd}{factor. indicates a zero R&D stock.}
  \item{randd}{numeric. R&D stock (in millions of dollars).}
  \item{employment}{numeric. employment (in 1000s).}
  \item{sp500}{factor. membership of firms in the S&P500 index.}
  \item{tobinq}{numeric. Tobin's q.}
  \item{value}{numeric. stock market value.}
  \item{institutions}{numeric. proportion of stock owned by institutions.}
  \item{industry}{factor. four-digit industry code.}
  \item{year}{factor. estimation period.}
  \item{top1}{numeric. share of the largest institution.}
  \item{quasiindexed}{numeric. share of "quasi-indexed" institutional owners.}
  \item{nonquasiindexed}{numeric. share of "non-quasi-indexed"
    institutional owners.}
  \item{transient}{numeric. share of "transient" institutional owners.}
  \item{dedicated}{numeric. share of "dedicated" institutional owners.}
  \item{competition4}{numeric. FIXME.}
  \item{subsample}{factor. subsample for the replication of columns
    1--5 from Table 4 in Aghion et al. (2013).}
}
}

\details{
Aghion et al. (2013) combine several firm level panel datasets (e.g.,
USPTO, SEC and Compustat) to examine the role of institutional investors
in the governance of innovation.

The classification of institutions into "quasiindexed", "transient" and
"dedicated" follows Bushee (1998).
}

\source{
Data and online appendix of Aghion et al. (2013).
}

\references{
Aghion P, Van Reenen J, Zingales L (2013).
  \dQuote{Innovation and Institutional Ownership.}
  \emph{The American Economic Review}, \bold{103}(1), 277--304.
  \doi{10.1257/aer.103.1.277}
  
Berger S, Stocker H, Zeileis A (2016).
  \dQuote{Innovation and Institutional Ownership Revisited: An Empirical Investigation with Count Data Models.}
  \emph{Empirical Economics}. Forthcoming.
  \doi{10.1007/s00181-016-1118-0}

  Blundell R, Griffith R, Van Reenen J (1999).
  \dQuote{Market Share, Market Value and Innovation in a Panel of
    British Manufacturing Firms.}
  \emph{Review of Economic Studies}, {66}(3), 529--554.
  
Bushee B (1998).
  \dQuote{The Influence of Institutional Investors on Myopic R&D Investment Behavior.}
  \emph{Accounting Review}, \bold{73}(3), 655--679.
}

\examples{
## Poisson models from Table I in Aghion et al. (2013)
## load data set
data(InstInnovation)

## log-scale variable
InstInnovation$lograndd <- log(InstInnovation$randd)
InstInnovation$lograndd[InstInnovation$lograndd == -Inf] <- 0

## count regression formulas
f1 <- cites ~ institutions + log(capital/employment) + log(sales) + industry + year
f2 <- cites ~ institutions + log(capital/employment) + log(sales) +
industry + year + lograndd + drandd
f3 <- cites ~ institutions + log(capital/employment) + log(sales) +
industry + year + lograndd + drandd + dprecites + log(precites)

## Poisson models
tab_I_3_pois <- glm(f1, data = InstInnovation, family = poisson)
tab_I_4_pois <- glm(f2, data = InstInnovation, family = poisson)
tab_I_5_pois <- glm(f3, data = InstInnovation, family = poisson)

## one-way clustered covariances
(vCL_I_3 <- vcovCL(tab_I_3_pois, cluster = InstInnovation$company))
(vCL_I_4 <- vcovCL(tab_I_4_pois, cluster = InstInnovation$company))
(vCL_I_5 <- vcovCL(tab_I_5_pois, cluster = InstInnovation$company))

## replication of columns 3 to 5 from Table I in Aghion et al. (2013)
lmtest::coeftest(tab_I_3_pois, vcov = vCL_I_3)
lmtest::coeftest(tab_I_4_pois, vcov = vCL_I_4)
lmtest::coeftest(tab_I_5_pois, vcov = vCL_I_5)

## model innovation as a dichotomous variable
InstInnovation$cites_binary <- ifelse(InstInnovation$cites == 0, 0, 1)

## logistic regression formulas
f4 <- update(f1, cites_binary ~ .)
f5 <- update(f2, cites_binary ~ .)
f6 <- update(f3, cites_binary ~ .)

## logistic regression
tab_I_3_logit <- glm(f4, data = InstInnovation, family = binomial)
tab_I_4_logit <- glm(f5, data = InstInnovation, family = binomial)
tab_I_5_logit <- glm(f6, data = InstInnovation, family = binomial)

## one-way clustered covariances
(vCL_I_3l <- vcovCL(tab_I_3_logit, cluster = InstInnovation$company))
(vCL_I_4l <- vcovCL(tab_I_4_logit, cluster = InstInnovation$company))
(vCL_I_5l <- vcovCL(tab_I_5_logit, cluster = InstInnovation$company))

##
lmtest::coeftest(tab_I_3_logit, vcov = vCL_I_3l)
lmtest::coeftest(tab_I_4_logit, vcov = vCL_I_4l)
lmtest::coeftest(tab_I_5_logit, vcov = vCL_I_5l)

}

\keyword{datasets}
