\name{Investment}
\alias{Investment}
\title{US Investment Data}
\description{
US data for fitting an investment equation.
}
\usage{data(Investment)}
\format{
An annual time series from 1963 to 1982 with 7 variables.
 \describe{
    \item{GNP}{nominal gross national product (in billion USD),}
    \item{Investment}{nominal gross private domestic investment (in billion USD),}
    \item{Price}{price index, implicit price deflator for GNP,}
    \item{Interest}{interest rate, average yearly discount rate
      charged by the New York Federal Reserve Bank,}
    \item{RealGNP}{real GNP (= GNP/Price),}
    \item{RealInv}{real investment (= Investment/Price),}
    \item{RealInt}{approximation to the real interest rate
      (= Interest - 100 * diff(Price)/Price).}
 }
}

\source{Table 15.1 in Greene (1993)}

\references{
Greene W.H. (1993), \emph{Econometric Analysis}, 2nd edition.
  Macmillan Publishing Company, New York.

Executive Office of the President (1984), \emph{Economic Report of the
  President}. US Government Printing Office, Washington, DC.
}

\examples{
## Willam H. Greene, Econometric Analysis, 2nd Ed.
## Chapter 15
## load data set, p. 411, Table 15.1
data(Investment)

## fit linear model, p. 412, Table 15.2
fm <- lm(RealInv ~ RealInt + RealGNP, data = Investment)
summary(fm)

## visualize residuals
plot(ts(residuals(fm), start = 1964),
  type = "b", pch = 19, ylim = c(-35, 35), ylab = "Residuals")
sigma <- sqrt(sum(residuals(fm)^2)/fm$df.residual) ## maybe used df = 26 instead of 16 ??
abline(h = c(-2, 0, 2) * sigma, lty = 2)

if(require(lmtest)) {
## Newey-West covariances, Example 15.3
coeftest(fm, vcov = NeweyWest(fm, lag = 4))
## Note, that the following is equivalent:
coeftest(fm, vcov = kernHAC(fm, kernel = "Bartlett", bw = 5, prewhite = FALSE, adjust = FALSE))

## Durbin-Watson test, p. 424, Example 15.4
dwtest(fm, data = Investment)

## Breusch-Godfrey test, p. 427, Example 15.6
bgtest(fm, order = 4, data = Investment)
}
}

\keyword{datasets}
