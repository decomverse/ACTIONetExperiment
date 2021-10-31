#ifndef RCPP_UTIL_H
#define RCPP_UTIL_H

namespace ACTIONetExperiment
{
  template <typename T>
  // Rcpp::NumericVector arma2vec(const T &x);
  Rcpp::NumericVector fast_row_sums(SEXP &A);
  Rcpp::NumericVector fast_column_sums(SEXP &A);
  Rcpp::NumericVector fast_row_max(SEXP &A);
  arma::sp_mat bind_mats_sparse(SEXP &X1, SEXP &X2, int dim);
  arma::mat bind_mats_dense(SEXP &X1, SEXP &X2, int dim);

  template <typename T>
  T bind_mats_generic(const SEXP &X1, const SEXP &X2, int dim);
  //
  template <class T1, class T2>
  bool kv_pair_less(const std::pair<T1, T2> &x, const std::pair<T1, T2> &y);

  void csr_sort_indices_inplace(Rcpp::IntegerVector &Ap, Rcpp::IntegerVector &Aj, Rcpp::NumericVector &Ax);
  void csc_sort_indices_inplace(Rcpp::IntegerVector &Ap, Rcpp::IntegerVector &Ai, Rcpp::NumericVector &Ax);
}
#endif
