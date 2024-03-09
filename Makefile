# Common variables
CC = mpicxx
CFLAGS = -std=c++11 -Wall -O0 -pedantic
LIBS = -lboost_program_options -lblas -llapack -lscalapack-openmpi

# Source file directory
SRC_DIR = src
# Object file directory
BUILD_DIR = build
# Executable to be created
OUTPUT = solver
# List of source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
# List of object files
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Number of processors
NP ?= 1

all: $(BUILD_DIR) $(OUTPUT)

$(OUTPUT): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

run: $(OUTPUT)
	mpiexec -np $(NP) ./$(OUTPUT)

clean:
	rm -rf $(BUILD_DIR)/*.o $(OUTPUT)

.PHONY: all clean