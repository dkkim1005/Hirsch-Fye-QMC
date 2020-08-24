#pragma once
#include <exception>
#include <string>

extern "C"
{
  void dger_(const int * m, const int * n, const double * alpha, const double * x, const int * incx,
    const double * y, const int * incy, double * a, const int * lda);

  void dgesv_(const int * n, const int * nrhs, double * a, const int * lda, int * ipiv, double * b,
    const int * ldb, int * info);
}

namespace blas
{
  inline void ger(const int m, const int n, const double alpha, const double * x, const double * y, double * a)
  {
    const int inc = 1;
    dger_(&m, &n, &alpha, x, &inc, y, &inc, a, &m);
  }
} // end namespace blas

namespace lapack
{
  inline void gesv(const int n, const int nrhs, double * a, double * b)
  {
    std::vector<int> ipiv(n);
    int info;
    dgesv_(&n, &nrhs, a, &n, ipiv.data(), b, &n, &info);
    if (info != 0)
      throw std::runtime_error("# dgesv info: " + std::to_string(info));
  }
} // end namespace lapack
