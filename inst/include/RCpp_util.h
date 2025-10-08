#ifndef RCPP_UTIL_H
#define RCPP_UTIL_H

namespace ACTIONetExperiment
{
  template <class T1, class T2>
  bool kv_pair_less(const std::pair<T1, T2> &x, const std::pair<T1, T2> &y);

  void csr_sort_indices_inplace(Rcpp::IntegerVector &Ap, Rcpp::IntegerVector &Aj, Rcpp::NumericVector &Ax);
  void csc_sort_indices_inplace(Rcpp::IntegerVector &Ap, Rcpp::IntegerVector &Ai, Rcpp::NumericVector &Ax);
}
#endif
