\name{cdfDensityBounded}
\alias{cdfDensityBounded}
\alias{quantileDensityBounded}

\title{
Cumulative distribution and quantiles of univariate model-based mixture density estimation for bounded data
}

\description{
Compute the cumulative density function (cdf) or quantiles of a one-dimensional density for bounded data estimated via transformation-based approach for Gaussian mixtures using \code{\link{densityMclustBounded}}.}

\usage{
cdfDensityBounded(object, data, ngrid = 100, \dots)
quantileDensityBounded(object, p, \dots)
}

\arguments{
  \item{object}{a \code{densityMclustBounded} model object.}
  \item{data}{a numeric vector of evaluation points.}
  \item{ngrid}{the number of points in a regular grid to be used as evaluation points if no \code{data} are provided.}
  \item{p}{a numeric vector of probabilities.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{The cdf is evaluated at points given by the optional argument \code{data}. If not provided, a regular grid of length \code{ngrid} for the evaluation points is used. 

The quantiles are computed using bisection linear search algorithm.
}

\value{
\code{cdfDensityBounded} returns a list of \code{x} and \code{y} values providing, respectively, the evaluation points and the estimated cdf. 

\code{quantileDensityBounded} returns a vector of quantiles.
}

\author{Luca Scrucca}

\seealso{
  \code{\link{densityMclustBounded}}, 
  \code{\link{plot.densityMclustBounded}}.
}

\examples{
\donttest{
# univariate case with lower bound
x <- rchisq(200, 3)
dens <- densityMclustBounded(x, lbound = 0)

xgrid <- seq(-2, max(x), length=1000)
cdf <- cdfDensityBounded(dens, xgrid)
str(cdf)
plot(xgrid, pchisq(xgrid, df = 3), type = "l", xlab = "x", ylab = "CDF")
lines(cdf, col = 4, lwd = 2)

q <- quantileDensityBounded(dens, p = c(0.01, 0.1, 0.5, 0.9, 0.99))
cbind(quantile = q, cdf = cdfDensityBounded(dens, q)$y)
plot(cdf, type = "l", col = 4, xlab = "x", ylab = "CDF")
points(q, cdfDensityBounded(dens, q)$y, pch = 19, col = 4)

# univariate case with lower & upper bounds
x <- rbeta(200, 5, 1.5)
dens <- densityMclustBounded(x, lbound = 0, ubound = 1)

xgrid <- seq(-0.1, 1.1, length=1000)
cdf <- cdfDensityBounded(dens, xgrid)
str(cdf)
plot(xgrid, pbeta(xgrid, 5, 1.5), type = "l", xlab = "x", ylab = "CDF")
lines(cdf, col = 4, lwd = 2)

q <- quantileDensityBounded(dens, p = c(0.01, 0.1, 0.5, 0.9, 0.99))
cbind(quantile = q, cdf = cdfDensityBounded(dens, q)$y)
plot(cdf, type = "l", col = 4, xlab = "x", ylab = "CDF")
points(q, cdfDensityBounded(dens, q)$y, pch = 19, col = 4)
}
}
