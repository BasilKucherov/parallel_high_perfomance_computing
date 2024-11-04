#include "matrix_generator.hpp"

#include <iomanip>
#include <sstream>

namespace {
enum class GridCellType { NORMAL, SPLITTED_LEFT, SPLITTED_RIGHT };

inline int GetSplittedCount(int gridIndex, int K1, int K2) {
  // Number of splitted cells BEFORE this index (not including this index)
  return gridIndex / (K1 + K2) * K2 +
         ((gridIndex % (K1 + K2)) > K1 ? (gridIndex % (K1 + K2)) - K1 : 0);
}

inline int GetGridCellIndex(int index, int K1, int K2) {
  return index - (index / (K1 + 2 * K2) * K2 +
                  (((index % (K1 + 2 * K2)) > K1)
                       ? (((index % (K1 + 2 * K2)) - K1 + 1) / 2)
                       : 0));
}

inline bool IsGridCellSplitted(int gridIndex, int K1, int K2) {
  return (gridIndex % (K1 + K2)) >= K1;
}

inline GridCellType GetCellType(int index, int K1, int K2) {
  auto indexInGroup = index % (K1 + 2 * K2);
  if (indexInGroup < K1) {
    return GridCellType::NORMAL;
  }

  if ((indexInGroup - K1) % 2 == 0) {
    return GridCellType::SPLITTED_LEFT;
  }

  return GridCellType::SPLITTED_RIGHT;
}

inline int GetRealCellIndex(int index, int K1, int K2) {
  return index + GetSplittedCount(index, K1, K2);
}

inline int GetBlockConnectionsCount(int Nx, int Ny) {
  return Nx * Ny * 2 - (Nx + Ny);
}

inline int GetJAIndex(int index, int i, int j, int K1, int K2, int Nx, int Ny) {
  int leftBlockConnectionsCount =
      ((i + 1) * j != 0) ? (2 * GetBlockConnectionsCount(i + 1, j) + i + 1 +
                            (((i + 1) != Ny) ? j : 0))
                         : 0;
  int topBlockConnectionsCount =
      (i * (Nx - j) != 0) ? (2 * GetBlockConnectionsCount(i, Nx - j) + Nx - j +
                             ((j != 0) ? i : 0))
                          : 0;

  auto totalConnectionsCount = leftBlockConnectionsCount +
                               topBlockConnectionsCount +
                               GetSplittedCount(i * Nx + j, K1, K2) * 2 + index;

  const auto cellType = GetCellType(index, K1, K2);
  if (cellType == GridCellType::SPLITTED_RIGHT) {
    totalConnectionsCount += (i != 0);
    totalConnectionsCount += (j != 0);
    totalConnectionsCount++;
  }
  return totalConnectionsCount;
}
} // namespace

namespace NMathUtils {
std::string GetMatrixPortrait(const int *IA, const int *JA, int N) {
  std::stringstream ss;
  for (int i = 0; i < N; ++i) {
    ss << std::setw(2) << i << ":*";
    auto jIndex = IA[i];
    auto maxJIndex = IA[i + 1] - 1;
    for (int j = 0; j < N; ++j) {
      if (JA[jIndex] == j && jIndex <= maxJIndex) {
        ss << "X ";
        jIndex++;
      } else {
        ss << "  ";
      }
    }
    ss << "*\n";
  }
  return ss.str();
}

void GenerateMatrixCSR(int Nx, int Ny, int K1, int K2, int *&IA, int *&JA,
                       int &N) {
  int initialCellCount = Nx * Ny;
  int splitCellCount = GetSplittedCount(initialCellCount, K1, K2);

  N = initialCellCount + splitCellCount;
  int connectionsCount = 2 * (2 * Nx * Ny - (Nx + Ny) + splitCellCount) + N;

  IA = new int[N + 1]();
  JA = new int[connectionsCount];

#pragma omp parallel for
  for (int cellIndex = 0; cellIndex < N; ++cellIndex) {

    const auto gridCellType = GetCellType(cellIndex, K1, K2);
    const int gridCellIndex = GetGridCellIndex(cellIndex, K1, K2);

    int i = gridCellIndex / Nx;
    int j = gridCellIndex % Nx;

    auto currentJAIndex = GetJAIndex(cellIndex, i, j, K1, K2, Nx, Ny);
    IA[cellIndex] = currentJAIndex;

    if (i != 0 && gridCellType != GridCellType::SPLITTED_RIGHT) {
      const auto topGridCellIndex = gridCellIndex - Nx;
      auto topCellIndex = GetRealCellIndex(topGridCellIndex, K1, K2) +
                          IsGridCellSplitted(topGridCellIndex, K1, K2);

      JA[currentJAIndex++] = topCellIndex;
    }

    if (j != 0 || gridCellType == GridCellType::SPLITTED_RIGHT) {
      JA[currentJAIndex++] = cellIndex - 1;
    }

    JA[currentJAIndex++] = cellIndex;

    if (j != Nx - 1 || gridCellType == GridCellType::SPLITTED_LEFT) {
      JA[currentJAIndex++] = cellIndex + 1;
    }

    if (i != Ny - 1 && gridCellType != GridCellType::SPLITTED_LEFT) {
      const auto bottomGridCellIndex = gridCellIndex + Nx;
      auto bottomCellIndex = GetRealCellIndex(bottomGridCellIndex, K1, K2);
      JA[currentJAIndex++] = bottomCellIndex;
    }
  }

  IA[N] = connectionsCount;
}
} // namespace NMathUtils
