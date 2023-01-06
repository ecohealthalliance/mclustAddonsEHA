\name{logsumexp}
\alias{logsumexp}
\alias{softmax}

\title{
Compute log-sum-exp and softmax functions.
}

\description{
Efficient implementation (via Rcpp) of log-sum-exp and softmax functions. 
}

\usage{
logsumexp(x, v = NULL)
softmax(x, v = NULL)
}

\arguments{
  \item{x}{a matrix of dimension \eqn{(n \times k)}{(n x k)}.}
  \item{v}{an optional vector of length \eqn{k}{k}. If not provided a vector 
           of zeros is used.}
}

\details{Given the matrix \eqn{x}{x} \code{logsumexp()} calculates for each row 
\eqn{x_{[i,]}}{x[i,]} the log-sum-exp function computed as

\deqn{m_{[i]} + \log\sum(\exp(x_{[i,]} + v) - m_{[i]})}{m[i] + log(sum(exp(x[i,] + v) - m[i]))}

where \eqn{m_{[i]} = \max(x_{[i,]} + v)}{m[i] = max(x[i,] + v) }.

\code{softmax()} calculates for each row \eqn{x_{[i,]}}{x[i,]} the softmax (aka 
multinomial logistic) function 

\deqn{\frac{\exp( x_{[i,c]} + v_{[c]} )}{\sum_{j=1}^k\exp( x_{[i,j]} + v_{[j]})}}{exp( x_[i,c] + v_[c]) / (sum_{j=1}^k exp( x_[i,j] + v_[j]) )}
}

\value{\code{logsumexp()} returns a vector of values of length equal to the number of rows of \eqn{x}{x}.
\code{softmax()} returns a matrix of values of the same dimension as \eqn{x}{x}.}

\author{Luca Scrucca}

\examples{
x = matrix(rnorm(15), 5, 3)
v = log(c(0.5, 0.3, 0.2))
logsumexp(x, v)
(z = softmax(x, v))
rowSums(z)
}
