\name{densityMclustBounded.diagnostic}
\alias{densityMclustBounded.diagnostic}

\title{Diagnostic plots for \code{mclustDensityBounded} estimation}

\description{
Diagnostic plots for density estimation of bounded data via transformation-based approach of Gaussian mixtures.
Only available for the one-dimensional case.
}

\usage{
densityMclustBounded.diagnostic(object, type = c("cdf", "qq"), 
                                col = c("black", "black"), 
                                lwd = c(2,1), lty = c(1,1), 
                                legend = TRUE, grid = TRUE, 
                                \dots)
}

\arguments{
  \item{object}{An object of class \code{'mclustDensityBounded'} obtained from a call to \code{\link{densityMclustBounded}} function.}

  \item{type}{The type of graph requested:
    \describe{
      \item{\code{"cdf"} =}{a plot of the estimated CDF versus the empirical distribution function.}
      \item{\code{"qq"} =}{a Q-Q plot of sample quantiles versus the quantiles obtained from the inverse of the estimated cdf.}
    }
  }
  
  \item{col}{A pair of values for the color to be used for plotting, respectively, the estimated CDF and the empirical cdf.}

  \item{lwd}{A pair of values for the line width to be used for plotting, respectively, the estimated CDF and the empirical cdf.}

  \item{lty}{A pair of values for the line type to be used for plotting, respectively, the estimated CDF and the empirical cdf.}

  \item{legend}{A logical indicating if a legend must be added to the plot of fitted CDF vs the empirical CDF.}

  \item{grid}{A logical indicating if a \code{\link{grid}} should be added to the plot.}

  \item{\dots}{Additional arguments.}
}

\details{
The two diagnostic plots for density estimation in the one-dimensional case are discussed in Loader (1999, pp- 87-90).
}

\value{No return value, called for side effects.}

\references{
  Loader C. (1999), Local Regression and Likelihood. New York, Springer.
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
plot(dens, x, what = "diagnostic")
# or
densityMclustBounded.diagnostic(dens, type = "cdf")
densityMclustBounded.diagnostic(dens, type = "qq")

# univariate case with lower & upper bounds
x <- rbeta(200, 5, 1.5)
dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
plot(dens, x, what = "diagnostic")
# or
densityMclustBounded.diagnostic(dens, type = "cdf")
densityMclustBounded.diagnostic(dens, type = "qq")
}
}
