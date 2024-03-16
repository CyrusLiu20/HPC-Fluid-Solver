# Build system for 2D Lid Driven Cavity Incompressible Flow
# all: 			builds the main solver executable
# unittests: 	builds the test cases executable for the involved classes
# run_test: 	executes the test executable
# run: 			executes the main solver executable
# doc:			generates documentation using Doxygen
# clean:		cleans object files, executables, test report, documentation, and results

# Common variables
CC = mpicxx
MPI_CC = mpiexec
CFLAGS = -std=c++11 -Wall -O0 -pedantic
LIBS = -lboost_program_options -lblas -llapack -lscalapack-openmpi 
TEST_LIBS = -lboost_unit_test_framework


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

# Test report file to be generated
TEST_REPORT = test_report.txt

# Doxygen file generated
DOCX = Doxyfile

# parameters to run the main solver executable (can be overwritten)
Np ?= 1 # Number of processors
Lx ?= 1.0 # Length of the domain in the x-direction
Ly ?= 1.0 # Length of the domain in the y-direction
dt ?= 2e-4 # Time step
T  ?= 1.0 # Total simulation time
Re ?= 10 # Reynolds number
Nx ?= 9 # Number of grid points in x-direction
Ny ?= 9 # Number of grid points in y-direction
verbose ?= false # Verbosity


# Builds the main solver executable
all: $(BUILD_DIR) $(OUTPUT) doc

$(OUTPUT): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Builds the test cases executable for the involved classes
unittests: $(TEST_BUILD_DIR) $(TEST_OUTPUT)

$(TEST_BUILD_DIR):
	mkdir -p $(TEST_BUILD_DIR)

$(TEST_BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_OUTPUT): $(TEST_OBJS) $(TEST_SRC_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) $(TEST_LIBS)

# Executes the test executable
run_test: $(TEST_OUTPUT)
	$(MPI_CC) -np $(Np) ./$(TEST_OUTPUT) --log_level=test_suite --report_level=detailed --output_format=HRF --log_sink=$(TEST_REPORT)

# Executes the main solver executable
run: $(OUTPUT)
	$(MPI_CC) -np $(Np) ./$(OUTPUT) --Lx $(Lx) --Ly $(Ly) --dt $(dt) --T $(T) --Re $(Re) --Nx $(Nx) --Ny $(Ny) --verbose $(verbose)

doc:
	doxygen -g $(DOCX)

clean:
	rm -rf $(BUILD_DIR)/*.o $(OUTPUT) $(TEST_BUILD_DIR)/*.o $(TEST_OUTPUT) $(DOCX) $(TEST_REPORT) $(RES_DIR)/*.txt

.PHONY: all clean