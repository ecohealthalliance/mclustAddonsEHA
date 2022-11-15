#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector logsumexp_Rcpp(NumericMatrix x, NumericVector a)
{
// Efficiently computes log-sum-exp(x+a)
// x = matrix (n x d)
// a = vector (d)

  double n = x.nrow();
  double d = x.ncol();
  NumericVector lse(n);
  NumericVector xa(d);
  double m;
  // 
  for (int i = 0; i < n; i++)
  {
    xa = x.row(i) + a;
    m = max(xa);
    lse[i] = m + log(sum(exp(xa-m)));
  }
  return lse;
}

/***

x = structure(c(-4.19768936334846, -23.3334911845962, -5.36851445858848, 
-15.2896460085004, -3.06018772423303, -13.2857737610833, -3.69968442181734, 
-5.33468420156765, -22.954092839643, -3.03420360101199, -23.3405056884397, 
-2.6395810621981, -16.4338853853632, -3.23305725595493, -50.7400647615373, 
-7.76487486677727, -57.9522847161203, -26.7659048640944, -2.41249310267583, 
-44.9591733534474), dim = c(10L, 2L))
pro = c(0.644072589572232, 0.355927410427768)

microbenchmark::microbenchmark(
  logsumexp_Rcpp(x, log(pro)),
  logsumexp_Rcpp(sweep(x, MARGIN = 2, STATS = log(pro), FUN = "+"), rep(0,2)),
  apply(sweep(x, MARGIN = 2, STATS = log(pro), FUN = "+"), 1, mclust:::logsumexp))
*/


// [[Rcpp::export]]
NumericMatrix softmax_Rcpp(NumericMatrix x, NumericVector a)
{
// Efficiently computes softmax function based on log-sum-exp(x+a)
// x = matrix (n x d)
// a = vector (d)

  double n = x.nrow();
  double d = x.ncol();
  NumericVector lse = logsumexp_Rcpp(x, a);
  NumericVector xa(d);
  NumericMatrix z(n,d);
  // 
  for (int i = 0; i < n; i++)
  {
    xa = x.row(i) + a;
    z(i, _) = exp(xa- lse(i));
  }
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

microbenchmark::microbenchmark(
  softmax_Rcpp(x, log(pro)),
  {
    xx = sweep(x, MARGIN = 2, STATS = log(pro), FUN = "+")
    exp(sweep(xx, MARGIN = 1, FUN = "-", STATS = apply(xx, 1, mclust:::logsumexp)))
  })

*/
