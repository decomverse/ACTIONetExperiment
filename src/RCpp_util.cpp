
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

vec roll_var(vec &X)
{
  const uword n_max = X.n_elem;
  double xbar = 0, M = 0;
  vec out(n_max);
  double *x = X.begin(), *o = out.begin();

  for (uword n = 1; n <= n_max; ++n, ++x, ++o)
  {
    double tmp = (*x - xbar);
    xbar += (*x - xbar) / n;
    M += tmp * (*x - xbar);
    if (n > 1L)
      *o = M / (n - 1.);
  }

  if (n_max > 0)
    out[0] = std::numeric_limits<double>::quiet_NaN();

  return out;
}

Rcpp::NumericVector fast_row_sums(SEXP &A)
{
  vec sum_vec;
  if (Rf_isS4(A))
  {
    sp_mat X = as<arma::sp_mat>(A);
    sum_vec = zeros(X.n_rows);

    sp_mat::const_iterator it = X.begin();
    sp_mat::const_iterator it_end = X.end();
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

Rcpp::NumericVector fast_column_sums(SEXP &A)
{
  vec sum_vec;
  if (Rf_isS4(A))
  {
    sp_mat X = as<arma::sp_mat>(A);
    sum_vec = zeros(X.n_cols);

    sp_mat::const_iterator it = X.begin();
    sp_mat::const_iterator it_end = X.end();
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

Rcpp::NumericVector fast_row_max(SEXP &A)
{
  vec sum_vec;
  if (Rf_isS4(A))
  {
    sp_mat X = as<arma::sp_mat>(A);
    sum_vec = zeros(X.n_rows);

    sp_mat::const_iterator it = X.begin();
    sp_mat::const_iterator it_end = X.end();
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

// Adapted from
// https://github.com/GreenleafLab/MPAL-Single-Cell-2019/blob/master/scRNA_02_Cluster_Disease_w_Reference_v1.R
Rcpp::NumericVector computeSparseRowVariances(IntegerVector j,
                                              NumericVector val,
                                              NumericVector rm, int n)
{
  const int nv = j.size();
  const int nm = rm.size();
  Rcpp::NumericVector rv(nm);
  Rcpp::NumericVector rit(nm);
  int current;
  // Calculate RowVars Initial
  for (int i = 0; i < nv; ++i)
  {
    current = j(i) - 1;
    rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
    rit(current) = rit(current) + 1;
  }
  // Calculate Remainder Variance
  for (int i = 0; i < nm; ++i)
  {
    rv(i) = rv(i) + (n - rit(i)) * rm(i) * rm(i);
  }
  rv = rv / (n - 1);
  return (rv);
}

sp_mat bind_sparse_mats(sp_mat &A, sp_mat &B, int dim = 0)
{
  sp_mat C;
  if (dim == 0)
  {
    C = join_cols(A, B);
  }
  else if (dim == 1)
  {
    C = join_rows(A, B);
  }
  else
  {
    stderr_stop("Invalid dim");
  }
  return (C);
}

template <class T1, class T2>
bool kv_pair_less(const std::pair<T1, T2> &x, const std::pair<T1, T2> &y)
{
  return x.first < y.first;
}

// [[Rcpp::export]]
void csr_sort_indices_inplace(IntegerVector &Ap, IntegerVector &Aj,
                              NumericVector &Ax)
{
  int n_row = Ap.size() - 1;
  std::vector<std::pair<int, double> > temp;

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
void csc_sort_indices_inplace(IntegerVector &Ap, IntegerVector &Ai,
                              NumericVector &Ax)
{
  int n_col = Ap.size() - 1;

  std::vector<std::pair<int, double> > temp;
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
