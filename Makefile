# Common variables
CC = mpicxx
MPI_CC = mpiexec
CFLAGS = -std=c++11 -Wall -O0 -pedantic
LIBS = -lboost_program_options -lblas -llapack -lscalapack-openmpi

# Source file directory
SRC_DIR = src
# Object file directory
BUILD_DIR = build
# Results file directory
RES_DIR = results
# Executable to be created
OUTPUT = solver
# Doxygen file generated
DOCX = Doxyfile
# List of source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
# List of object files
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

NP ?= 1 # Number of processors
Lx ?= 1.0 # Length of the domain in the x-direction
Ly ?= 1.0 # Length of the domain in the y-direction
Nx ?= 9 # Number of grid points in x-direction
Ny ?= 9 # Number of grid points in y-direction

all: $(BUILD_DIR) $(OUTPUT) doc

$(OUTPUT): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

run: $(OUTPUT)
	$(MPI_CC) -np $(NP) ./$(OUTPUT) --Lx $(Lx) --Ly $(Ly) --Nx $(Nx) --Ny $(Ny)

doc:
	doxygen -g $(DOCX)

clean:
	rm -rf $(BUILD_DIR)/*.o $(OUTPUT) $(DOCX) $(RES_DIR)/*.txt

.PHONY: all clean