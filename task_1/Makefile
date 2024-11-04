CXX = clang++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -fopenmp -Iinclude -lgtest
GENERATES = run
TRASH = *.o

SRC_DIR = src
UTILS_DIR = $(SRC_DIR)/utils
BUILD_DIR = build
SRCS = $(SRC_DIR)/main.cpp $(UTILS_DIR)/arg_parser.cpp $(UTILS_DIR)/utils.cpp $(SRC_DIR)/matrix_generator.cpp $(SRC_DIR)/filler.cpp $(SRC_DIR)/solver.cpp

OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)

all: $(GENERATES)

$(BUILD_DIR)/%.o: %.cpps
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

run: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:
	rm -rf $(GENERATES) $(TRASH) $(BUILD_DIR)
