#include <vector>

namespace utils {
enum class ErrorCode {
  OK,
  INVALID_ARGUMENTS,
};

void PrintHelp(const char *programName);
} // namespace utils