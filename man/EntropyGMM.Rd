\name{EntropyGMM}
\alias{EntropyGMM}
\alias{EntropyGMM.densityMclust}
\alias{EntropyGMM.Mclust}
\alias{EntropyGMM.densityMclustBounded}
\alias{EntropyGMM.matrix}
\alias{EntropyGMM.data.frame}
\alias{EntropyGauss}
\alias{nats2bits}
\alias{bits2nats}

\title{Gaussian mixture-based estimation of entropy}

\description{
Compute an estimate of the (differential) entropy from a Gaussian Mixture Model (GMM) fitted using the \emph{mclust} package.
}

\usage{

EntropyGMM(object, \dots) 

\method{EntropyGMM}{densityMclust}(object, \dots)
\method{EntropyGMM}{Mclust}(object, \dots)
\method{EntropyGMM}{densityMclustBounded}(object, \dots)
\method{EntropyGMM}{matrix}(object, \dots)
\method{EntropyGMM}{data.frame}(object, \dots)

EntropyGauss(sigma)

nats2bits(x)
bits2nats(x)
}

\arguments{

\item{object}{An object of class \code{'Mclust'}, \code{'densityMclust'}, or \code{'densityMclustBounded'}, obtained by fitting a Gaussian mixture via, respectively, \code{\link[mclust]{Mclust}()}, \code{\link[mclust]{densityMclust}()}, and \code{\link{densityMclustBounded}()}. 
If a \code{'matrix'} or \code{'data.frame'} is provided as input, a GMM using the provided data is estimated preliminary to computing the entropy. In this case further arguments can be provided to control the fitted model (e.g. number of mixture components and/or covariances decomposition).}

\item{sigma}{A symmetric covariance matrix.}

\item{x}{A vector of values.}

\item{\dots}{Further arguments passed to or from other methods.}
}

\value{
\code{EntropyGMM()} returns an estimate of the entropy based on a estimated Gaussian mixture model (GMM) fitted using the \emph{mclust} package. If a matrix of data values is provided, a GMM is preliminary fitted to the data and then the entropy computed.

\code{EntropyGauss()} returns the entropy for a multivariate Gaussian distribution with covariance matrix \code{sigma}.

\code{nats2bits()} and \code{bits2nats()} convert values in nats to bits, and viceversa. Information-theoretic quantities have different units depending on the base of the logarithm used: nats are expressed in base-2 logarithms, whereas bits in natural logarithms.
}

\seealso{
\code{\link[mclust]{Mclust}}, \code{\link[mclust]{densityMclust}}.
}

\author{Luca Scrucca}

\references{
Robin S. and Scrucca L. (2023) Mixture-based estimation of entropy. \emph{Computational Statistics & Data Analysis}, 177, 107582. https://doi.org/10.1016/j.csda.2022.107582
}

\examples{
\donttest{
X = iris[,1:4]
mod = densityMclust(X, plot = FALSE)
h = EntropyGMM(mod)
h
bits2nats(h)
EntropyGMM(X)
}}

