
#include <RcppArmadillo.h>

using namespace arma;
using namespace std;

#include <RCpp_util.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG
#define DYNSCHED

#define SYS_THREADS_DEF (std::thread::hardware_concurrency() - 2)

#define stdout_printf Rprintf
#define stderr_printf REprintf
#define stderr_stop stop
#define FLUSH R_FlushConsole()

template <typename T>
Rcpp::NumericVector arma2vec(const T &x)
{
  return Rcpp::NumericVector(x.begin(), x.end());
}


// [[Rcpp::export]]
Rcpp::NumericVector fast_row_sums(SEXP &A)
{
  arma::vec sum_vec;
  if (Rf_isS4(A))
  {
    arma::sp_mat X = as<arma::sp_mat>(A);
    sum_vec = zeros(X.n_rows);

    arma::sp_mat::const_iterator it = X.begin();
    arma::sp_mat::const_iterator it_end = X.end();
    for (; it != it_end; ++it)
    {
      sum_vec[it.row()] += (*it);
    }
  }
  else
  {
    mat X = as<arma::mat>(A);
    sum_vec = sum(X, 1);
  }

  return (arma2vec(sum_vec));
}

// [[Rcpp::export]]
Rcpp::NumericVector fast_column_sums(SEXP &A)
{
  arma::vec sum_vec;
  if (Rf_isS4(A))
  {
    arma::sp_mat X = as<arma::sp_mat>(A);
    sum_vec = zeros(X.n_cols);

    arma::sp_mat::const_iterator it = X.begin();
    arma::sp_mat::const_iterator it_end = X.end();
    for (; it != it_end; ++it)
    {
      sum_vec[it.col()] += (*it);
    }
  }
  else
  {
    mat X = as<arma::mat>(A);
    sum_vec = trans(sum(X, 0));
  }

  return (arma2vec(sum_vec));
}

// [[Rcpp::export]]
Rcpp::NumericVector fast_row_max(SEXP &A)
{
  arma::vec sum_vec;
  if (Rf_isS4(A))
  {
    arma::sp_mat X = as<arma::sp_mat>(A);
    sum_vec = zeros(X.n_rows);

    arma::sp_mat::const_iterator it = X.begin();
    arma::sp_mat::const_iterator it_end = X.end();
    for (; it != it_end; ++it)
    {
      sum_vec[it.row()] = std::max(sum_vec[it.row()], (*it));
    }
  }
  else
  {
    mat X = as<arma::mat>(A);
    sum_vec = max(X, 1);
  }

  return (arma2vec(sum_vec));
}

template <typename T>
T bind_mats_generic(const T &A, const T &B, int dim)
{
  T C;

  if (dim == 0)
  {
    C = join_vert(A, B);
  }
  else if (dim == 1)
  {
    C = join_horiz(A, B);
  }
  else
  {
    stderr_stop("Invalid dim");
  }
  return (C);
}

template arma::mat bind_mats_generic<arma::mat>(const arma::mat &A, const arma::mat &B, int dim);
template arma::sp_mat bind_mats_generic<arma::sp_mat>(const arma::sp_mat &A, const arma::sp_mat &B, int dim);

// [[Rcpp::export]]
arma::sp_mat bind_mats_sparse(SEXP &X1, SEXP &X2, int dim = 0){
  arma::sp_mat X = bind_mats_generic<arma::sp_mat>(as<arma::sp_mat>(X1), as<arma::sp_mat>(X2), dim);
  return (X);
}
// [[Rcpp::export]]
arma::mat bind_mats_dense(SEXP &X1, SEXP &X2, int dim = 0){
  arma::mat X = bind_mats_generic<arma::mat>(as<arma::mat>(X1), as<arma::mat>(X2), dim);
  return (X);
}

template <class T1, class T2>
bool kv_pair_less(const std::pair<T1, T2> &x, const std::pair<T1, T2> &y)
{
  return x.first < y.first;
}

// [[Rcpp::export]]
void csr_sort_indices_inplace(Rcpp::IntegerVector &Ap, Rcpp::IntegerVector &Aj,
                              Rcpp::NumericVector &Ax)
{
  int n_row = Ap.size() - 1;
  std::vector<std::pair<int, double>> temp;

  for (int i = 0; i < n_row; i++)
  {
    int row_start = (int)Ap[i];
    int row_end = (int)Ap[i + 1];
    int len = row_end - row_start;

    temp.resize(len);
    bool is_sorted = true;
    for (int jj = row_start, n = 0; jj < row_end; jj++, n++)
    {
      temp[n].first = (int)Aj(jj);
      temp[n].second = Ax(jj);
      if ((jj < (row_end - 1)) && (Aj(jj + 1) < Aj(jj)))
      {
        is_sorted = false;
      }
    }
    if (is_sorted)
      continue;

    std::sort(temp.begin(), temp.begin() + len, kv_pair_less<int, double>);
    for (int jj = row_start, n = 0; jj < row_end; jj++, n++)
    {
      Aj(jj) = temp[n].first;
      Ax(jj) = temp[n].second;
    }
  }
}

// [[Rcpp::export]]
void csc_sort_indices_inplace(Rcpp::IntegerVector &Ap, Rcpp::IntegerVector &Ai,
                              Rcpp::NumericVector &Ax)
{
  int n_col = Ap.size() - 1;

  std::vector<std::pair<int, double>> temp;
  for (int i = 0; i < n_col; i++)
  {
    int col_start = (int)Ap[i];
    int col_end = (int)Ap[i + 1];
    int len = col_end - col_start;

    temp.resize(len);
    bool is_sorted = true;
    for (int jj = col_start, n = 0; jj < col_end; jj++, n++)
    {
      temp[n].first = (int)Ai(jj);
      temp[n].second = Ax(jj);
      if ((jj < (col_end - 1)) && (Ai(jj + 1) < Ai(jj)))
      {
        is_sorted = false;
      }
    }
    if (is_sorted)
      continue;

    std::sort(temp.begin(), temp.begin() + len, kv_pair_less<int, double>);
    for (int jj = col_start, n = 0; jj < col_end; jj++, n++)
    {
      Ai(jj) = temp[n].first;
      Ax(jj) = temp[n].second;
    }
  }
}
