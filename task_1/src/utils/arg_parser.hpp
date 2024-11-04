#include <cstddef>

namespace utils {
struct Arguments {
  bool Help;
  bool DebugOutput;
  int Threads;
  int Nx;
  int Ny;
  int K1;
  int K2;
};

class ArgParser {
public:
  ArgParser(int argc, char *argv[]);
  Arguments Parse();

private:
  int Argc_;
  char **Argv_;
};
} // namespace utils
