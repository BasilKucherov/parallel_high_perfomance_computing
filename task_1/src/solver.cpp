#include "solver.hpp"

#include <cmath>
#include <cstring>
#include <iostream>

namespace {
void AdamarsMult(int N, const double *X, const double *Y, double *Res) {
#pragma omp for
  for (int i = 0; i < N; ++i) {
    Res[i] = X[i] * Y[i];
  }
}

void SpMV(int N, const int *IA, const int *JA, const double *A, const double *X,
          double *Res) {
#pragma omp for
  for (int i = 0; i < N; ++i) {
    Res[i] = 0;
    for (int j = IA[i]; j < IA[i + 1]; ++j) {
      Res[i] += A[j] * X[JA[j]];
    }
  }
}

void DotProduct(int N, const double *X, const double *Y, double &Res) {
  double local_res = 0.0;

#pragma omp for
  for (int i = 0; i < N; ++i) {
    local_res += X[i] * Y[i];
  }

#pragma omp atomic
  Res += local_res;
}

void Axpy(int N, double a, const double *X, const double *Y, double *Res) {
#pragma omp for
  for (int i = 0; i < N; ++i) {
    Res[i] = a * X[i] + Y[i];
  }
}
} // namespace
namespace MathUtils {
void Solve(int N, const int *IA, const int *JA, const double *A,
           const double *B, double Eps, int MaxIter, double *&X, int &Iter,
           double &Res, bool Debug) {
  X = new double[N];

  double *AInvDiag = new double[N];
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    auto jIndex = IA[i];
    while (jIndex < IA[i + 1] && JA[jIndex] != i) {
      ++jIndex;
    }
    AInvDiag[i] = 1.0 / A[jIndex];
  }

  double *r = new double[N];
  double *p = new double[N];
  double *q = new double[N];
  double *z = new double[N];
  double rho, rhoPrev;
  double alpha;

  double tmp;
  double *tmpVec = new double[N];

  std::memcpy(r, B, N * sizeof(double));

  Iter = 0;
#pragma omp parallel
  do {
#pragma omp barrier
    AdamarsMult(N, r, AInvDiag, z);
    rho = 0;
#pragma omp barrier
    DotProduct(N, r, z, rho);

    if (Iter == 0) {
      std::memcpy(p, z, N * sizeof(double));
    } else {
#pragma omp barrier
      double beta = rho / rhoPrev;
#pragma omp barrier
      Axpy(N, beta, p, z, p);
    }
#pragma omp barrier
    SpMV(N, IA, JA, A, p, q);

    alpha = 0;
#pragma omp barrier
    DotProduct(N, p, q, alpha);
#pragma omp barrier

    alpha = rho / alpha;
#pragma omp barrier
    Axpy(N, alpha, p, X, X);
#pragma omp barrier

    Axpy(N, -alpha, q, r, r);

    if (Debug) {
      SpMV(N, IA, JA, A, X, tmpVec);
      Axpy(N, -1, tmpVec, B, tmpVec);
      DotProduct(N, tmpVec, tmpVec, tmp);
      Res = std::sqrt(tmp);

      std::cout << "Iter: " << Iter << std::endl;
      std::cout << "\tResidual: " << Res << std::endl;
      std::cout << "\trho: " << rho << std::endl;
      std::cout << "\trho > eps^2: " << (rho > Eps * Eps) << std::endl;
      std::cout << std::endl;
    }

    rhoPrev = rho;
  } while (rho > Eps * Eps && Iter++ < MaxIter);

  delete[] tmpVec;
  delete[] AInvDiag;
  delete[] r;
  delete[] p;
  delete[] q;
  delete[] z;
}

} // namespace MathUtils
