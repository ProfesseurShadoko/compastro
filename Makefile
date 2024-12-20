# Compiler and flags
CXX := g++
CXXFLAGS := -Wall -O2 -I/usr/include/eigen3 -MMD -MP -DCL_TARGET_OPENCL_VERSION=300

# Library dependencies
LIB := -lfftw3 -lm -lsfml-graphics -lsfml-window -lsfml-system -lOpenCL

# Project structure
SRC_DIR := src
BUILD_DIR := build

# Collect all .cpp files in the src directory
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
DEPFILES := $(OBJECTS:.o=.d)

# The target executable
TARGET := $(BUILD_DIR)/output.exe

# Default target: compile and run
all: $(TARGET) run

# Link the object files to create the executable
$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILD_DIR)
	@echo "Linking object files into $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LIB)

# Rule to compile source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	@echo "Compiling $< into $@..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Include the dependency files
-include $(DEPFILES)


# Run the program
run: $(TARGET)
	@echo "Running the program...\n"
	./$(TARGET)

# Clean build directory
clean:
	@echo "Cleaning up..."
	rm -rf $(BUILD_DIR)
