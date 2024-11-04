#include <vector>

namespace utils {
enum class ErrorCode {
  OK,
  INVALID_ARGUMENTS,
};

void PrintHelp(const char *programName);

struct MatrixPortrait {
  ;
  int *JA;
  int N;
};

void GenerateMatrixCSR(int Nx, int Ny, int K1, int K2, int *&IA, int *&JA,
                       int &N);

std::string GetMatrixPortrait(const int *IA, const int *JA, int N);

} // namespace utils