#include "solver.hpp"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <omp.h>

namespace MathUtils {
void AdamarsMult(int N, const double *X, const double *Y, double *Res) {
#pragma omp for
  for (int i = 0; i < N; ++i) {
    Res[i] = X[i] * Y[i];
  }
}

void SpMV(int N, const int *IA, const int *JA, const double *A, const double *X,
          double *Res) {
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    Res[i] = 0;
    for (int j = IA[i]; j < IA[i + 1]; ++j) {
      Res[i] += A[j] * X[JA[j]];
    }
  }
}

void DotProduct(int N, const double *X, const double *Y, double &Res) {
  Res = 0;
#pragma omp parallel for reduction(+ : Res)
  for (int i = 0; i < N; ++i) {
    Res += X[i] * Y[i];
  }
}

void Axpy(int N, double a, const double *X, const double *Y, double *Res) {
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    Res[i] = a * X[i] + Y[i];
  }
}

void Solve(int N, const int *IA, const int *JA, const double *A,
           const double *B, double Eps, int MaxIter, double *&X, int &Iter,
           double &Res, bool Debug) {
  double sumSpMVTime = 0;
  int spMVCount = 0;
  double sumDotProductTime = 0;
  int dotProductCount = 0;
  double sumAxpyTime = 0;
  int axpyCount = 0;
  double startTime;

  double *AInvDiag = new double[N];
#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    auto jIndex = IA[i];
    while (jIndex < IA[i + 1] && JA[jIndex] != i) {
      ++jIndex;
    }
    AInvDiag[i] = 1.0 / A[jIndex];
  }

  X = new double[N];
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
  do {
    AdamarsMult(N, r, AInvDiag, z);

    ++dotProductCount;
    startTime = omp_get_wtime();
    DotProduct(N, r, z, rho);
    sumDotProductTime += omp_get_wtime() - startTime;

    if (Iter == 0) {
      std::memcpy(p, z, N * sizeof(double));
    } else {
      double beta = rho / rhoPrev;
      ++axpyCount;
      startTime = omp_get_wtime();
      Axpy(N, beta, p, z, p);
      sumAxpyTime += omp_get_wtime() - startTime;
    }
    ++spMVCount;
    startTime = omp_get_wtime();
    SpMV(N, IA, JA, A, p, q);
    sumSpMVTime += omp_get_wtime() - startTime;

    ++dotProductCount;
    startTime = omp_get_wtime();
    DotProduct(N, p, q, alpha);
    sumDotProductTime += omp_get_wtime() - startTime;

    alpha = rho / alpha;
    ++axpyCount;
    startTime = omp_get_wtime();
    Axpy(N, alpha, p, X, X);
    sumAxpyTime += omp_get_wtime() - startTime;

    ++axpyCount;
    startTime = omp_get_wtime();
    Axpy(N, -alpha, q, r, r);
    sumAxpyTime += omp_get_wtime() - startTime;

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

  std::cout << std::setprecision(3)
            << "AVG/TOTAL SpMV:       " << sumSpMVTime / spMVCount << " / "
            << sumSpMVTime << std::endl;
  std::cout << std::setprecision(3)
            << "AVG/TOTAL DotProduct: " << sumDotProductTime / dotProductCount
            << " / " << sumDotProductTime << std::endl;
  std::cout << std::setprecision(3)
            << "AVG/TOTAL Axpy:       " << sumAxpyTime / axpyCount << " / "
            << sumAxpyTime << std::endl;

  delete[] tmpVec;
  delete[] AInvDiag;
  delete[] r;
  delete[] p;
  delete[] q;
  delete[] z;
}

} // namespace MathUtils
