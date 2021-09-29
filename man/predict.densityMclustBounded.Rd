\name{predict.densityMclustBounded}
\alias{predict.densityMclustBounded}

\title{Model-based mixture density estimation for bounded data}

\description{
Compute density estimation for univariate and multivariate bounded data based on Gaussian finite mixture models estimated by \code{\link{densityMclustBounded}}.
}

\usage{
\method{predict}{densityMclustBounded}(object, newdata, 
        what = c("dens", "cdens", "z"), 
        logarithm = FALSE, \dots)
}

\arguments{
\item{object}{An object of class \code{"densityMclustBounded"} resulting from a call to \code{\link{densityMclustBounded}}.}

\item{newdata}{A numeric vector, matrix, or data frame of observations. If missing the density is computed for the input data obtained from the call to \code{\link{densityMclustBounded}}.}

\item{what}{A character string specifying what to retrieve: \code{"dens"} returns a vector of values for the mixture density; \code{"cdens"} returns a matrix of component densities for each mixture component (along the columns); \code{"z"} returns a matrix of component posterior probabilities.}

\item{logarithm}{A logical value indicating whether or not the logarithm of the densities/probabilities should be returned.}

\item{\dots}{Further arguments passed to or from other methods.}
}

\value{
Returns a vector or a matrix of values evaluated at \code{newdata} depending on the argument \code{what} (see above).
}

\seealso{
\code{\link{densityMclustBounded}}, 
\code{\link{plot.densityMclustBounded}}.
}

\author{Luca Scrucca}

\references{
Scrucca L. (2019) A transformation-based approach to Gaussian mixture density estimation for bounded data. \emph{Biometrical Journal}, 61:4, 873–888. https://doi.org/10.1002/bimj.201800174
}

\examples{\donttest{
y <- sample(0:1, size = 200, replace = TRUE, prob = c(0.6, 0.4))
x <- y*rchisq(200, 3) + (1-y)*rchisq(200, 10)
dens <- densityMclustBounded(x, lbound = 0)
summary(dens)
plot(dens, what = "density", data = x, breaks = 11)

xgrid <- seq(0, max(x), length = 201)
densx <- predict(dens, newdata = xgrid, what = "dens")
cdensx <- predict(dens, newdata = xgrid, what = "cdens")
cdensx <- sweep(cdensx, MARGIN = 2, FUN = "*", dens$parameters$pro)
plot(xgrid, densx, type = "l", lwd = 2)
matplot(xgrid, cdensx, type = "l", col = 3:4, lty = 2:3, lwd = 2, add = TRUE)

z <- predict(dens, newdata = xgrid, what = "z")
matplot(xgrid, z, col = 3:4, lty = 2:3, lwd = 2, ylab = "Posterior probabilities")
}}
