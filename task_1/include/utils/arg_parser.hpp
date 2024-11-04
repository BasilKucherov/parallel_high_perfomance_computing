#include <cstddef>

namespace utils {
struct Arguments {
  int Nx;
  int Ny;
  int K1;
  int K2;
  int Threads;
  bool Help;
  bool DebugOutput;
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
