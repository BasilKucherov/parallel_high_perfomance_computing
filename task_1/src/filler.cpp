#include "filler.hpp"

#include <cmath>

namespace {
const auto PI = std::atan(1) * 4;
constexpr auto COEFF = 1.984;
} // namespace

namespace NMathUtils {
void FillMatrixCSR(int N, const int *IA, const int *JA, double *&A,
                   double *&B) {
  A = new double[IA[N]]();
  B = new double[N];

#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    double sum = 0;
    int diagIndex = 0;

    for (int j = IA[i]; j < IA[i + 1]; ++j) {
      if (JA[j] != i) {
        const auto value = std::cos(i * JA[j] + i + JA[j]);
        A[j] = value;
        sum += std::abs(value);
      } else {
        diagIndex = j;
      }
    }

    A[diagIndex] = COEFF * sum;
  }

#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    B[i] = std::sin(i);
  }
}
} // namespace NMathUtils
