#include "utils/arg_parser.hpp"

#include <cstring>
#include <stdexcept>
#include <string>

namespace utils {
ArgParser::ArgParser(int argc, char *argv[]) : Argc_(argc), Argv_(argv) {}

Arguments ArgParser::Parse() {
  auto args = Arguments{0, 0, 0, 0, 1, false, false};
  size_t PositionalArgsProcessedCount = 0;

  for (int i = 1; i < Argc_; ++i) {
    if (Argv_[i][0] == '-') {
      if (std::strcmp(Argv_[i], "-h") == 0 ||
          std::strcmp(Argv_[i], "--help") == 0) {
        args.Help = true;
      } else if (std::strcmp(Argv_[i], "-d") == 0 ||
                 std::strcmp(Argv_[i], "--debug") == 0) {
        args.DebugOutput = true;
      } else if (std::strcmp(Argv_[i], "-t") == 0 ||
                 std::strcmp(Argv_[i], "--threads") == 0) {
        args.Threads = std::strtol(Argv_[++i], nullptr, 10);
      } else {
        throw std::invalid_argument("Unknown flag: " + std::string(Argv_[i]));
      }
    } else {
      if (PositionalArgsProcessedCount == 0) {
        args.Nx = std::strtol(Argv_[i], nullptr, 10);
      } else if (PositionalArgsProcessedCount == 1) {
        args.Ny = std::strtol(Argv_[i], nullptr, 10);
      } else if (PositionalArgsProcessedCount == 2) {
        args.K1 = std::strtol(Argv_[i], nullptr, 10);
      } else if (PositionalArgsProcessedCount == 3) {
        args.K2 = std::strtol(Argv_[i], nullptr, 10);
      } else {
        throw std::invalid_argument("Too many positional arguments");
      }

      ++PositionalArgsProcessedCount;
    }
  }

  if (!args.Help && PositionalArgsProcessedCount != 4) {
    throw std::invalid_argument("Invalid number of positional arguments");
  }

  return args;
}
} // namespace utils
