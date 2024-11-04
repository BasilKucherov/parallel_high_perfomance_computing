#include "utils/utils.hpp"

#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>

namespace utils {
void PrintHelp(const char *programName) {
  std::cout << "Usage: " << programName
            << " [-d | --debug] [-h | --help] [-t | --threads] Nx Ny K1 K2\n";
}
} // namespace utils