CXX = clang++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -g -fopenmp
GENERATES = run
TRASH = *.o

SRC_DIR = src
UTILS_DIR = $(SRC_DIR)/utils
BUILD_DIR = build
SRCS = $(SRC_DIR)/main.cpp $(UTILS_DIR)/arg_parser.cpp $(UTILS_DIR)/utils.cpp

OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)

all: $(GENERATES)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

run: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

clean:
	rm -rf $(GENERATES) $(TRASH) $(BUILD_DIR)
