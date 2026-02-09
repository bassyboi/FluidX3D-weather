CC = g++
CFLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic
LDFLAGS = -lm

SRC_DIR = .
PHYSICS_DIR = physics
IO_DIR = io
BUILD_DIR = build
BIN_DIR = bin

# Source files
SRCS = main.cpp config.cpp grid.cpp terrain.cpp \
       $(PHYSICS_DIR)/solver.cpp $(PHYSICS_DIR)/thermal.cpp $(PHYSICS_DIR)/cloud.cpp \
       $(IO_DIR)/gfs.cpp $(IO_DIR)/radiosonde.cpp $(IO_DIR)/aws.cpp $(IO_DIR)/output.cpp

# Object files
OBJS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Test sources
TEST_CONFIG_SRC = tests/test_config.cpp config.cpp
TEST_SIM_SRC = tests/test_simulation.cpp config.cpp grid.cpp terrain.cpp \
               $(PHYSICS_DIR)/solver.cpp $(PHYSICS_DIR)/thermal.cpp $(PHYSICS_DIR)/cloud.cpp \
               $(IO_DIR)/gfs.cpp $(IO_DIR)/radiosonde.cpp $(IO_DIR)/aws.cpp $(IO_DIR)/output.cpp

TARGET = $(BIN_DIR)/weather

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(OBJS) -o $@ $(CFLAGS) $(LDFLAGS)
	@echo "Build successful: $@"

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CC) -c $< -o $@ $(CFLAGS)

# Tests
test: test_config test_simulation
	@echo "All tests passed!"

test_config: $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_CONFIG_SRC))
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $(BIN_DIR)/$@ $(CFLAGS) $(LDFLAGS)
	cd .. && ./weather/$(BIN_DIR)/$@

test_simulation: $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SIM_SRC))
	@mkdir -p $(BIN_DIR)
	$(CC) $^ -o $(BIN_DIR)/$@ $(CFLAGS) $(LDFLAGS)
	cd .. && ./weather/$(BIN_DIR)/$@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) output/

run: $(TARGET)
	cd .. && ./weather/$(TARGET) weather/config.toml
