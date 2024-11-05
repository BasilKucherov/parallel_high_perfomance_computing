#include "test_solver.hpp"
#include "solver.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <omp.h>

namespace {
const double EPS = 1e-6;
const int MIN_THREADS = 1;
const int MAX_THREADS = 8;
} // namespace

namespace TestUtils {
bool AllClose(double *A, double *B, int N, double eps = 1e-6) {
  for (int i = 0; i < N; ++i) {
    if (std::abs(A[i] - B[i]) > eps) {
      return false;
    }
  }
  return true;
}

void MatrixToCSR(int N, const double *M, int *&IA, int *&JA, double *&A) {
  int NonZero = 0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (M[i * N + j] != 0) {
        ++NonZero;
      }
    }
  }

  IA = new int[N + 1];
  JA = new int[NonZero];
  A = new double[NonZero];

  int k = 0;
  for (int i = 0; i < N; ++i) {
    IA[i] = k;
    for (int j = 0; j < N; ++j) {
      if (M[i * N + j] != 0) {
        A[k] = M[i * N + j];
        JA[k++] = j;
      }
    }
  }
  IA[N] = NonZero;
}

void MV(int N, const double *M, const double *X, double *&Res) {
  Res = new double[N];
  for (int i = 0; i < N; ++i) {
    Res[i] = 0;
    for (int j = 0; j < N; ++j) {
      Res[i] += M[i * N + j] * X[j];
    }
  }
}
} // namespace TestUtils

namespace Test {
void TestAdamarsMult() {
  std::cout << "Testing AdamarsMult:" << std::endl;
  int N = 10;
  double *X = new double[N];
  double *Y = new double[N];
  double *Res = new double[N];
  double *Expected = new double[N];

  for (int i = 0; i < N; ++i) {
    X[i] = std::sin(i);
    Y[i] = std::cos(i);
    Expected[i] = X[i] * Y[i];
  }

  for (int threads = MIN_THREADS; threads <= MAX_THREADS; ++threads) {
    std::cout << "\tTesting with " << threads << " threads ";
    omp_set_num_threads(threads);
    MathUtils::AdamarsMult(N, X, Y, Res);
    assert(TestUtils::AllClose(Res, Expected, N, EPS));
    std::cout << "OK" << std::endl;
  }

  delete[] X;
  delete[] Y;
  delete[] Res;
  delete[] Expected;
}

void TestSpMV() {
  std::cout << "Testing SpMV:" << std::endl;
  int N = 10;
  double *M = new double[N * N];
  double *X = new double[N];
  double *Res = new double[N];

  // Fill matrix
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      M[i * N + j] = std::sin(i) * std::cos(j);
    }
  }

  int *IA;
  int *JA;
  double *A;
  TestUtils::MatrixToCSR(N, M, IA, JA, A);

  double *Expected;
  TestUtils::MV(N, M, X, Expected);

  for (int threads = MIN_THREADS; threads <= MAX_THREADS; ++threads) {
    std::cout << "\tTesting with " << threads << " threads ";
    omp_set_num_threads(threads);
    MathUtils::SpMV(N, IA, JA, A, X, Res);
    assert(TestUtils::AllClose(Res, Expected, N, EPS));
    std::cout << "OK" << std::endl;
  }

  delete[] M;
  delete[] X;
  delete[] Res;
  delete[] Expected;
}

void TestDotProduct() {
  std::cout << "Testing DotProduct:" << std::endl;
  int N = 10;
  double *X = new double[N];
  double *Y = new double[N];
  double Expected = 0;
  double Res = 0;

  for (int i = 0; i < N; ++i) {
    X[i] = std::sin(i);
    Y[i] = std::cos(i);
    Expected += X[i] * Y[i];
  }

  for (int threads = MIN_THREADS; threads <= MAX_THREADS; ++threads) {
    Res = 0;
    std::cout << "\tTesting with " << threads << " threads ";
    omp_set_num_threads(threads);
    MathUtils::DotProduct(N, X, Y, Res);
    assert(std::abs(Res - Expected) < EPS);
    std::cout << "OK" << std::endl;
  }
}
} // namespace Test
