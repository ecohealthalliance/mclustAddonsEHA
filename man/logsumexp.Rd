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
logsumexp(x, a = NULL)
softmax(x, a = NULL)
}

\arguments{
  \item{x}{a matrix of dimension \eqn{(n \times k)}{(n x k)}.}
  \item{a}{an optional vector of length \eqn{k}{k}. If not provided a vector of 
           zeros is used.}
}

\details{Given the matrix \eqn{x}{x} \code{logsumexp()} calculates for each row 
\eqn{x_{[i,]}}{x[i,]} the log-sum-exp function computed as

\deqn{m_{[i]} + \log\sum(\exp(x_{[i,]} + a) - m_{[i]})}{m[i] + log(sum(exp(x[i,] + a) - m[i]))}

where \eqn{m_{[i]} = \max(x_{[i,]} + a)}{m[i] = max(x[i,] + a) }.

\code{softmax()} calculates for each row \eqn{x_{[i,]}}{x[i,]} the softmax (aka 
multinomial logistic) function 

\deqn{\dfrac{\exp( x_{[i,c]} + a_{[c]} )}{\sum_{j=1}^k\exp( x_{[i,j]} + a_{[j]})}}{exp( x_[i,c] + a_[c]) / (sum_{j=1}^k exp( x_[i,j] + a_[j]) )}
}

\author{Luca Scrucca}

\examples{
x = matrix(rnorm(15), 5, 3)
a = log(c(0.5, 0.3, 0.2))
logsumexp(x, a)
(z = softmax(x, a))
rowSums(z)
}
