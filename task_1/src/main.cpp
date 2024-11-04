#include "matrix_generator.hpp"
#include "utils/arg_parser.hpp"
#include "utils/utils.hpp"

#include "omp.h"
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <optional>

using utils::ErrorCode;

int main(int argc, char *argv[]) {
  auto argParser = utils::ArgParser(argc, argv);
  utils::Arguments args;

  try {
    args = argParser.Parse();
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    utils::PrintHelp(argv[0]);
    return static_cast<int>(ErrorCode::INVALID_ARGUMENTS);
  }

  if (args.Help) {
    utils::PrintHelp(argv[0]);
    return static_cast<int>(ErrorCode::OK);
  }

  if (args.Nx <= 0 || args.Ny <= 0 || args.K1 < 0 || args.K2 < 0) {
    std::cerr << "Error: Nx, Ny must be greater than 0, K1, K2 must be grater "
                 "or equal to 0"
              << std::endl;
    return static_cast<int>(ErrorCode::INVALID_ARGUMENTS);
  }

  if (args.Threads <= 0) {
    std::cerr << "Error: Threads must be greater than 0" << std::endl;
    return static_cast<int>(ErrorCode::INVALID_ARGUMENTS);
  } else {
    omp_set_num_threads(args.Threads);
  }

  int N, *IA, *JA;

  auto start = std::chrono::high_resolution_clock::now();
  CSR::GenerateMatrixCSR(args.Nx, args.Ny, args.K1, args.K2, IA, JA, N);
  auto finish = std::chrono::high_resolution_clock::now();

  std::cout << "Time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(finish -
                                                                     start)
                   .count()
            << "ms" << std::endl;

  std::cout << "N: " << N << std::endl;
  std::cout << "IA[N]: " << IA[N] << std::endl;

  if (args.DebugOutput) {
    const auto portrait = CSR::GetMatrixPortrait(IA, JA, N);
    std::cout << "\n" << portrait << std::endl;
  }

  return static_cast<int>(ErrorCode::OK);
}
