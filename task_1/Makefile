CXX = clang++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -fopenmp -Iinclude
TEST_CXXFLAGS = $(CXXFLAGS) -Itests
GENERATES = run run_test
TRASH = *.o

SRC_DIR = src
UTILS_DIR = $(SRC_DIR)/utils
BUILD_DIR = build
TEST_DIR = tests

MAIN_SRCS = $(SRC_DIR)/main.cpp
SRCS = $(UTILS_DIR)/arg_parser.cpp $(UTILS_DIR)/utils.cpp $(SRC_DIR)/matrix_generator.cpp $(SRC_DIR)/filler.cpp $(SRC_DIR)/solver.cpp
MAIN_TEST_SRCS = $(TEST_DIR)/main.cpp
TEST_SRCS = $(TEST_DIR)/test_solver.cpp

OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
TEST_OBJS = $(TEST_SRCS:%.cpp=$(BUILD_DIR)/%.o)
all: $(GENERATES)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

run: $(OBJS) $(MAIN_SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $(MAIN_SRCS) $(OBJS)

run_test: $(TEST_OBJS) $(MAIN_TEST_SRCS)
	$(CXX) $(TEST_CXXFLAGS) -o $@ $(MAIN_TEST_SRCS) $(TEST_OBJS) $(OBJS)

clean:
	rm -rf $(GENERATES) $(TRASH) $(BUILD_DIR)
