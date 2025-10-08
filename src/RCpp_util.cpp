
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
