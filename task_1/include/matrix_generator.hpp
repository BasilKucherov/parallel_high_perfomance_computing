#include <string>

namespace NMathUtils {
void GenerateMatrixCSR(int Nx, int Ny, int K1, int K2, int *&IA, int *&JA,
                       int &N);

std::string GetMatrixPortrait(const int *IA, const int *JA, int N);
} // namespace NMathUtils
