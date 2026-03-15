// scNPM_types.h
#ifndef SCNPM_TYPES_H
#define SCNPM_TYPES_H

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#define ARMA_NO_DEBUG
#define RCPP_USE_GLOBAL_STRINGSTREAM

#include <omp.h>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace arma;

inline vec rcpp_normpdf(const vec& x) {
  vec res(x.n_elem);
  for (uword i = 0; i < x.n_elem; i++) {
    res(i) = R::dnorm(x(i), 0.0, 1.0, 0);
  }
  return res;
}

inline vec rcpp_normcdf(const vec& x) {
  vec res(x.n_elem);
  for (uword i = 0; i < x.n_elem; i++) {
    res(i) = R::pnorm(x(i), 0.0, 1.0, 1, 0);
  }
  return res;
}

inline double rcpp_normpdf(double x) {
  return R::dnorm(x, 0.0, 1.0, 0);
}

inline double rcpp_normcdf(double x) {
  return R::pnorm(x, 0.0, 1.0, 1, 0);
}

#endif // SCNPM_TYPES_H
