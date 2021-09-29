#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//
// Range-Transformation
// 
// [[Rcpp::export]]
Rcpp::NumericVector 
  rangeTransform(NumericVector x, 
                 double lbound = NA_REAL,
                 double ubound = NA_REAL)
{
  NumericVector z = clone(x);
  
  if(is_finite(lbound) && is_finite(ubound))
  {
    z = (x - lbound)/(ubound - x);
  }
  else if(is_finite(lbound))
  {
    z = (x - lbound);
  }
  else if(is_finite(ubound)) 
  { 
    stop("upper bound only range transformation not available!");
  }

  return z;
}

/*** R
# set.seed(1)
# x = rchisq(5,3)
# rangeTransform_R(x)
# rangeTransform(x)
# rangeTransform_R(x, lbound = -1)
# rangeTransform(x, lbound = -1)
# rangeTransform_R(x, lbound = 0, ubound = 10)
# rangeTransform(x, lbound = 0, ubound = 10)
# 
# set.seed(123)
# x = rchisq(1e8,3)
# system.time(rangeTransform_R(x, lbound = 0, ubound = 10))
# system.time(rangeTransform(x, lbound = 0, ubound = 10))
*/

//
// Power Box-Cox Transformation
// 
// [[Rcpp::export]]
Rcpp::NumericVector 
  powerTransform(NumericVector x, 
                 double lambda = 1,
                 double eps = 1e-3)
{
  // NumericVector x = clone(data);
  // LogicalVector isna = is_na(x);
  // x = x[!isna];
  NumericVector z = clone(x);
  
  if(is_true(any(z < 0)))
    { stop("data values must be strictly positive."); }
  
  if(std::abs(lambda) <= eps)
    { z = log(x); } 
  else 
    { z = (pow(x, lambda) - 1)/lambda; }

  return(z);
}

/*** R
# x = rchisq(10, 3)
# powerTransform_R(x, lambda = 1/3)
# powerTransform(x, lambda = 1/3)
# powerTransform_R(x, lambda = 0)
# powerTransform(x, lambda = 0)
# 
# set.seed(123)
# x = rchisq(1e8,3)
# system.time(powerTransform_R(x, lambda = 1/3))
# system.time(powerTransform(x, lambda = 1/3))
# system.time(powerTransform_R(x, lambda = 0))
# system.time(powerTransform(x, lambda = 0))
*/

//
// Transformation derivatives
// 

// [[Rcpp::export]]
Rcpp::NumericVector 
  rangepowerTransformDeriv_lb(NumericVector x,
                              double lambda = 1,
                              double lbound = NA_REAL,
                              double eps    = NA_REAL)
{
  if(ISNAN(lbound))
    { stop("lbound missing!"); }
  if(ISNAN(eps)) // (!is_finite(eps))
    { stop("ubound missing!"); }
  // Rcout << "eps = " << eps << std::endl;
  
  NumericVector tx = rangeTransform(x, lbound); // = (x-lbound)
  NumericVector dx = pow(tx, lambda - 1);
  if(lambda < 1)
    { // linear interpolation near the boundary
      double b = (lambda-1)*pow(eps, lambda-2);
      double a = pow(eps, lambda-1) - b*(lbound+eps);
      dx = ifelse(x < (lbound+eps), a+b*x, dx);
    }
        
  return(dx);
}

/*** R
# set.seed(1)
# x = rchisq(200, 3)
# dx_R = rangepowerTransformDeriv_R(x, lambda = 0.3, lbound = 0, eps = 1)
# dx = rangepowerTransformDeriv_lb(x, lambda = 0.3, lbound = 0, eps = 1)
# plot(x, dx_R)
# points(x, dx, pch = 20, col = 2)
# 
# system.time(123)
# x = rchisq(1e7, 3)
# system.time(rangepowerTransformDeriv_R(x, lambda = 0.3, lbound = 0, eps = 1))
# system.time(rangepowerTransformDeriv_lb(x, lambda = 0.3, lbound = 0, eps = 1))
*/ 

// [[Rcpp::export]]
Rcpp::NumericVector 
  rangepowerTransformDeriv_lub(NumericVector x,
                               double lambda = 1,
                               double lbound = NA_REAL,
                               double ubound = NA_REAL,
                               double eps    = NA_REAL,
                               double tol    = 1e-3)
{
  if(ISNAN(lbound))
    { stop("lbound missing!"); }
  if(ISNAN(ubound))
    { stop("ubound missing!"); }
  if(ISNAN(eps)) // (!is_finite(eps))
    { stop("ubound missing!"); }
  
  double machineEps = sqrt(std::numeric_limits<double>::epsilon());
  NumericVector dx(x.size(), 1.0);
  double b, a;
  
  if(std::abs(lambda) < tol)
    { dx = 1/(x-lbound)+1/(ubound-x); 
      // linear interpolation near the lower bound
      b = 1/pow(ubound - lbound - eps, 2) - 1/pow(eps, 2);
      a = 1/eps + 1/(ubound - lbound - eps) - b*(lbound+eps);
      dx = ifelse(x < (lbound+eps), a+b*x, dx);  
      // linear interpolation near the upper bound
      b = 1/pow(eps, 2) - 1/pow(ubound - lbound - eps, 2);
      a = 1/(ubound - lbound - eps) + 1/eps - b*(ubound - eps);
      dx = ifelse(x > (ubound-eps), pmax(machineEps, a+b*x), dx);  
    }
  else
    { NumericVector tx = rangeTransform(x, lbound, ubound); 
      // = (x-lbound)/(ubound-x)
      dx = pow(tx, lambda - 1)*(ubound - lbound);
      dx = dx/((ubound-x)* (ubound-x));
      // linear interpolation near the lower bound
      b = ( pow(eps/(ubound - lbound - eps), lambda-2) * 
            (eps/(ubound - lbound - eps) +1) * (lambda-1) + 
            2*pow(eps/(ubound - lbound - eps), lambda-1) ) * 
          (ubound - lbound)/pow(ubound - lbound - eps, 3);
      a = pow(eps/(ubound - lbound - eps), lambda-1) * 
          (ubound - lbound)/pow(ubound - lbound - eps, 2) - 
          b*(lbound+eps);
      dx = ifelse(x < (lbound+eps), pmax(machineEps, a+b*x), dx);  
      // linear interpolation near the upper bound
      b = ( pow((ubound - lbound - eps)/eps, lambda-2) * 
            ((ubound - lbound - eps)/eps +1) * (lambda-1) + 
            2*pow((ubound - lbound - eps)/eps, lambda-1) ) * 
          (ubound - lbound)/pow(eps, 3);
      a = pow((ubound - lbound - eps)/eps, lambda-1) * 
          (ubound - lbound)/pow(eps, 2) - b*(ubound - eps);
      dx = ifelse(x > (ubound-eps), pmax(machineEps, a+b*x), dx);  
    }
  return(dx);
}

/*** R
# set.seed(1)
# x = c(0.001, rbeta(5, 1/2, 1/2), 0.999)
# rangepowerTransformDeriv_R(x, lambda = 1, lbound = 0, ubound = 1)
# rangepowerTransformDeriv_lub(x, lambda = 1, lbound = 0, ubound = 1)
# rangepowerTransformDeriv_R(x, lambda = 0, lbound = 0, ubound = 1)
# rangepowerTransformDeriv_lub(x, lambda = 0, lbound = 0, ubound = 1)
# rangepowerTransformDeriv_R(x, lambda = -1, lbound = 0, ubound = 1)
# rangepowerTransformDeriv_lub(x, lambda = -1, lbound = 0, ubound = 1)
#   
# set.seed(1)
# x = c(0.01, rbeta(200, 3, 1), 0.999)
# hist(x)
# dx_R = rangepowerTransformDeriv_R(x, lambda = 0, lbound = 0, eps = 1)
# dx = rangepowerTransformDeriv_lb(x, lambda = 0, lbound = 0, eps = 1)
# plot(x, dx_R)
# points(x, dx, pch = 20, col = 2)
# 
# x = rbeta(1e7, 1/2, 1/2)
# system.time(rangepowerTransformDeriv_R(x,   lambda = 0.3, lbound = 0, ubound = 1))
# system.time(rangepowerTransformDeriv_lub(x, lambda = 0.3, lbound = 0, ubound = 1))
# system.time(rangepowerTransformDeriv_R(x,   lambda = 0, lbound = 0, ubound = 1))
# system.time(rangepowerTransformDeriv_lub(x, lambda = 0, lbound = 0, ubound = 1))
# system.time(rangepowerTransformDeriv_R(x,   lambda = -1, lbound = 0, ubound = 1))
# system.time(rangepowerTransformDeriv_lub(x, lambda = -1, lbound = 0, ubound = 1))
*/ 

//
// Maximum by rows for numerical matrices
// (much faster than using apply(..., 1, max))
//

// [[Rcpp::export]]
NumericVector rowMax_sugar(NumericMatrix X)
{
  int nrows = X.nrow();
  NumericVector maxs(nrows);
  for(int i = 0; i < nrows; i++)
     { maxs[i] = Rcpp::max( X(i,_) ); }
  return(maxs);
}

// [[Rcpp::export]]
NumericVector rowMax(arma::mat X)
{
  arma::colvec rmax = arma::max(X, 1);
  return NumericVector(rmax.begin(), rmax.end());
}

// [[Rcpp::export]]
NumericVector rowSum(arma::mat X)
{
  arma::colvec rsum = arma::sum(X, 1);
  return NumericVector(rsum.begin(), rsum.end());
}

// [[Rcpp::export]]
NumericVector colMax(arma::mat X)
{
  arma::rowvec cmax = arma::max(X, 0);
  return NumericVector(cmax.begin(), cmax.end());
}

// [[Rcpp::export]]
NumericVector colSum(arma::mat X)
{
  arma::rowvec csum = arma::sum(X, 0);
  return NumericVector(csum.begin(), csum.end());
}

/*** R
# library(microbenchmark)
# X = matrix(rchisq(1e5,10),1e4,10)
# 
# all.equal(apply(X, 1, max), rowMax(X))
# all.equal(apply(X, 1, max), rowMax_sugar(X))
# microbenchmark(apply(X, 1, max), rowMax(X), rowMax_sugar(X))
# 
# all.equal(apply(X, 1, sum), rowSums(X))
# all.equal(apply(X, 1, sum), rowSum(X))
# microbenchmark(apply(X, 1, sum), rowSums(X), rowSum(X))
# 
# all.equal(apply(X, 2, max), colMax(X))
# microbenchmark(apply(X, 2, max), colMax(X))
#
# all.equal(apply(X, 2, sum), colSums(X))
# all.equal(apply(X, 2, sum), colSum(X))
# microbenchmark(apply(X, 2, sum), colSums(X), colSum(X))
*/
