#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector logsumexp_Rcpp(arma::mat& x, arma::rowvec& v)
{
// Efficiently computes log-sum-exp(x+v)
// x = matrix (n x k)
// v = vector (k)
  int          n = x.n_rows;
  int          k = x.n_cols;
  arma::vec    lse(n);
  arma::rowvec xv(k);
  double       m = 0.0;
  //  
  for (int i = 0; i < n; i++)
  {
    xv = x.row(i) + v;
    m = max(xv);
    lse[i] = m + log(sum(exp(xv-m)));
  }
  //  
  return NumericVector(Rcpp::wrap(lse));
}

/***
x = structure(c(-4.19768936334846, -23.3334911845962, -5.36851445858848, 
-15.2896460085004, -3.06018772423303, -13.2857737610833, -3.69968442181734, 
-5.33468420156765, -22.954092839643, -3.03420360101199, -23.3405056884397, 
-2.6395810621981, -16.4338853853632, -3.23305725595493, -50.7400647615373, 
-7.76487486677727, -57.9522847161203, -26.7659048640944, -2.41249310267583, 
-44.9591733534474), dim = c(10L, 2L))
pro = c(0.644072589572232, 0.355927410427768)
logsumexp_Rcpp(x, log(pro))
apply(sweep(x, MARGIN = 2, STATS = log(pro), FUN = "+"), 1, mclust:::logsumexp)
 
microbenchmark::microbenchmark(
  "Rcpp" = logsumexp_Rcpp(x, log(pro)),
  "Rcpp+R" = logsumexp_Rcpp(sweep(x, MARGIN = 2, STATS = log(pro), FUN = "+"), rep(0.0,2)),
  "R" = apply(sweep(x, MARGIN = 2, STATS = log(pro), FUN = "+"), 1, mclust:::logsumexp)
)
*/


// [[Rcpp::export]]
arma::mat softmax_Rcpp(arma::mat& x, arma::rowvec& v)
{
// Efficiently computes softmax function based on log-sum-exp(x+v)
// x = matrix (n x k)
// v = vector (k)
  int          n = x.n_rows;
  int          k = x.n_cols;
  arma::vec    lse = logsumexp_Rcpp(x, v);
  arma::rowvec xv(k);
  arma::mat    z(n,k);
  //
  for (int i = 0; i < n; i++)
  {
    xv = x.row(i) + v;
    z.row(i) = exp(xv - lse(i));
  }
  //
  return z;
}

/***
x = structure(c(-4.19768936334846, -23.3334911845962, -5.36851445858848, 
-15.2896460085004, -3.06018772423303, -13.2857737610833, -3.69968442181734, 
-5.33468420156765, -22.954092839643, -3.03420360101199, -23.3405056884397, 
-2.6395810621981, -16.4338853853632, -3.23305725595493, -50.7400647615373, 
-7.76487486677727, -57.9522847161203, -26.7659048640944, -2.41249310267583, 
-44.9591733534474), dim = c(10L, 2L))
pro = c(0.644072589572232, 0.355927410427768)
softmax_Rcpp(x, log(pro))
softmax_R <- function(x, v)
{ 
 x = sweep(x, MARGIN = 2, STATS = v, FUN = "+")
 exp(sweep(x, MARGIN = 1, FUN = "-", STATS = apply(x, 1, mclust:::logsumexp)))
}
softmax_R(x, log(pro))
 
microbenchmark::microbenchmark(
  "Rcpp" = softmax_Rcpp(x, log(pro)),
  "R" = softmax_R(x, log(pro)))
*/
