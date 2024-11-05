#include "test_solver.hpp"

int main() {
  Test::TestAdamarsMult();
  Test::TestSpMV();
  Test::TestDotProduct();
  return 0;
}
