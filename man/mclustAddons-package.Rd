\name{mclustAddons-package}
\alias{mclustAddons-package}
\alias{mclustAddons}
\docType{package}

\title{Addons for the \pkg{mclust} package}

\description{
Extend the functionality of the \pkg{mclust} package for Gaussian finite mixture modeling by including: 
\itemize{
  \item density estimation for data with bounded support (Scrucca, 2019)
  \item modal clustering using MEM algorithm for Gaussian mixtures (Scrucca, 2021)
  \item entropy estimation via Gaussian mixture modeling (Robin & Scrucca, 2023)
 }
}

\details{
For a quick introduction to \pkg{mclustAddons} see the vignette \href{../doc/mclustAddons.html}{A quick tour of mclustAddons}.

See also:
\itemize{
  \item \code{\link{densityMclustBounded}()} for density estimation of bounded data;
  \item \code{\link{MclustMEM}()} for modal clustering;
  \item  \code{\link{EntropyGMM}()} for entropy estimation.
  }
}

\author{
Luca Scrucca.

Maintainer: Luca Scrucca \email{luca.scrucca@unipg.it}
}

\references{
Scrucca L. (2019) A transformation-based approach to Gaussian mixture density estimation for bounded data. \emph{Biometrical Journal}, 61:4, 873–888. https://doi.org/10.1002/bimj.201800174

Scrucca L. (2021) A fast and efficient Modal EM algorithm for Gaussian mixtures. \emph{Statistical Analysis and Data Mining}, 14:4, 305–314. https://doi.org/10.1002/sam.11527

Robin S. and Scrucca L. (2023) Mixture-based estimation of entropy. \emph{Computational Statistics & Data Analysis}, 177, 107582. https://doi.org/10.1016/j.csda.2022.107582
}

\keyword{package}
