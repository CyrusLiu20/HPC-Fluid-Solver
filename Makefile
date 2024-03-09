# Common variables
CC = mpicxx
MPI_CC = mpiexec
CFLAGS = -std=c++11 -Wall -O0 -pedantic
LIBS = -lboost_program_options -lblas -llapack -lscalapack-openmpi

# Source folder, Build folder, and Results directory
SRC_DIR = src
BUILD_DIR = build
RES_DIR = results
# Executable to be created
OUTPUT = solver
# List of source files and object files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Test source folder and Test build folder
TEST_DIR = test
TEST_BUILD_DIR = build/test
# Test executable
TEST_OUTPUT = testing
# List of test source files and test objective files
TESTS = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJS = $(patsubst $(TEST_DIR)/%.cpp,$(TEST_BUILD_DIR)/%.o,$(TESTS))
# Excluding the main object file when compiling the testing executable
EXCLUDE_SRCS = $(SRC_DIR)/LidDrivenCavitySolver.cpp
TEST_SRC_OBJS = $(filter-out $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(EXCLUDE_SRCS)),$(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS)))

# Doxygen file generated
DOCX = Doxyfile

NP ?= 1 # Number of processors
Lx ?= 1.0 # Length of the domain in the x-direction
Ly ?= 1.0 # Length of the domain in the y-direction
dt ?= 2e-4 # Time step
T  ?= 1.0 # Total simulation time
Re ?= 10 # Reynolds number
Nx ?= 9 # Number of grid points in x-direction
Ny ?= 9 # Number of grid points in y-direction


all: $(BUILD_DIR) $(OUTPUT) doc

$(OUTPUT): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)


test: $(TEST_BUILD_DIR) $(TEST_OUTPUT)

$(TEST_BUILD_DIR):
	mkdir -p $(TEST_BUILD_DIR)

$(TEST_BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_OUTPUT): $(TEST_OBJS) $(TEST_SRC_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)


run: $(OUTPUT)
	$(MPI_CC) -np $(NP) ./$(OUTPUT) --Lx $(Lx) --Ly $(Ly) --dt $(dt) --T $(T) --Re $(Re) --Nx $(Nx) --Ny $(Ny)

doc:
	doxygen -g $(DOCX)

clean:
	rm -rf $(BUILD_DIR)/*.o $(OUTPUT) $(TEST_BUILD_DIR)/*.o $(TEST_OUTPUT) $(DOCX) $(RES_DIR)/*.txt

.PHONY: all clean