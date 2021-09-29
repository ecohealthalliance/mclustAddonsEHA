\name{plot.densityMclustBounded}
\alias{plot.densityMclustBounded}

\title{Plotting method for model-based mixture density estimation for bounded data}

\description{
Plots for \code{mclustDensityBounded} objects.
}

\usage{
\method{plot}{densityMclustBounded}(x, what = c("BIC", "density", "diagnostic"), 
     data = NULL, \dots)
}

\arguments{
\item{x}{An object of class \code{"densityMclustBounded"} obtained from a call to \code{\link{densityMclustBounded}}.}

\item{what}{The type of graph requested:
\describe{
  \item{\code{"BIC"} =}{a plot of BIC values for the estimated models versus the number of components.}
  \item{\code{"density"} =}{a plot of estimated density; if \code{data} is also provided the density is plotted over data points.}
  \item{\code{"diagnostic"} =}{diagnostic plots (only available for the one-dimensional case).}
  }
}
  
\item{data}{Optional data points.}

\item{\dots}{Further available arguments.

\itemize{
\item For 1-dimensional data:

\code{hist.col = "lightgrey",  hist.border = "white", breaks = "Sturges"}

\item For 2-dimensional data:

\code{type = c("contour", "hdr", "image", "persp")} \cr
\code{transformation = c("none", "log", "sqrt")} \cr
\code{grid = 100, nlevels = 11, levels = NULL} \cr
\code{prob = c(0.25, 0.5, 0.75)} \cr
\code{col = grey(0.6), color.palette = blue2grey.colors} \cr
\code{points.col = 1, points.cex = 0.8, points.pch = 1}

\item For \eqn{d > 2}-dimensional data:

\code{type = c("contour", "hdr"), gap = 0.2} \cr
\code{grid = 100, nlevels = 11, levels = NULL} \cr
\code{prob = c(0.25, 0.5, 0.75)} \cr
\code{col = grey(0.6), color.palette = blue2grey.colors} \cr
\code{points.col = 1, points.cex = 0.8, points.pch = 1}
}
}
}

\value{No return value, called for side effects.}

\seealso{
\code{\link{densityMclustBounded}},
\code{\link{predict.densityMclustBounded}}.
}

\author{Luca Scrucca}

\references{
Scrucca L. (2019) A transformation-based approach to Gaussian mixture density estimation for bounded data. \emph{Biometrical Journal}, 61:4, 873–888. https://doi.org/10.1002/bimj.201800174
}

\examples{\donttest{
# univariate case with lower bound
x <- rchisq(200, 3)
dens <- densityMclustBounded(x, lbound = 0)
plot(dens, what = "BIC")
plot(dens, what = "density", data = x, breaks = 15)

# univariate case with lower & upper bound
x <- rbeta(200, 5, 1.5)
dens <- densityMclustBounded(x, lbound = 0, ubound = 1)
plot(dens, what = "BIC")
plot(dens, what = "density", data = x, breaks = 9)

# bivariate case with lower bounds
x1 <- rchisq(200, 3)
x2 <- 0.5*x1 + sqrt(1-0.5^2)*rchisq(200, 5)
x <- cbind(x1, x2)
dens <- densityMclustBounded(x, lbound = c(0,0))
plot(dens, what = "density")
plot(dens, what = "density", data = x)
plot(dens, what = "density", type = "hdr")
plot(dens, what = "density", type = "persp")
}}
