\name{GaussianMixtureMEM}
\alias{GaussianMixtureMEM}

\title{Modal EM algorithm for Gaussian Mixtures}

\description{
A function implementing a fast and efficient Modal EM algorithm for Gaussian mixtures.
}

\usage{
GaussianMixtureMEM(data, pro, mu, sigma,
                   control = list(eps = 1e-5, 
                                  maxiter = 1e3, 
                                  stepsize = function(t) 1-exp(-0.1*t),
                                  denoise = TRUE,
                                  alpha = 0.01,
                                  keep.path = FALSE),
                   \dots)
}

\arguments{

\item{data}{A numeric vector, matrix, or data frame of observations. Categorical variables are not allowed. If a matrix or data frame, rows correspond to observations (\eqn{n}) and columns correspond to variables (\eqn{d}).} 

\item{pro}{A \eqn{(G \times 1)}{(G x 1)} vector of mixing probabilities for a Gaussian mixture of \eqn{G} components.}

\item{mu}{A \eqn{(d \times G)}{(d x G)} matrix of component means for a \eqn{d}-variate Gaussian mixture of \eqn{G} components.}

\item{sigma}{A \eqn{(d \times d \times G)}{(d x d x G)} array of component covariance matrices for a \eqn{d}-variate Gaussian mixture of \eqn{G} components.}

\item{control}{A list of control parameters:

\describe{
\item{\code{eps, maxiter}}{numerical values setting the tolerance and the maximum number of iterations of the MEM algorithm;}
\item{\code{stepsize}}{a function controlling the step size of the MEM algorithm;}
\item{\code{denoise}}{a logical, if \code{TRUE} a denoising procedure is used when \eqn{d > 1} to discard all modes whose density is negligible;}
\item{\code{alpha}}{a numerical value used when \code{denoise = TRUE} for computing the hypervolume of central \eqn{(1-\alpha)100}{(1-alpha)100} region of a multivariate Gaussian:}
\item{\code{keep.path}}{a logical controlling whether or not the full paths to modes must be returned.}
}
}

\item{\dots}{Further arguments passed to or from other methods.}
}

\value{
Returns a list containing the following elements:
\item{n}{The number of input data points.}
\item{d}{The number of variables/features.}
\item{parameters}{The Gaussian mixture parameters.}
\item{iter}{The number of iterations of MEM algorithm.}
\item{nmodes}{The number of modes estimated by the MEM algorithm.}
\item{modes}{The coordinates of modes estimated by MEM algorithm.}
\item{path}{If requested, the coordinates of full paths to modes for each data point.}
\item{logdens}{The log-density at the estimated modes.}
\item{logvol}{The log-volume used for denoising (if requested).}
\item{classification}{The modal clustering classification of input data points.}
}

\seealso{
\code{\link{MclustMEM}}.
}

\author{Luca Scrucca}

\references{
Scrucca L. (2021) A fast and efficient Modal EM algorithm for Gaussian mixtures. \emph{Statistical Analysis and Data Mining}, 14:4, 305–314. https://doi.org/10.1002/sam.11527
}
