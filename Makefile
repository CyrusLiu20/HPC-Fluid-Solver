## Common variables
CC = mpicxx
FLAGS = -std=c++11 -Wall -O0 -pedantic
SOURCEDIR = src/
LIBS = -lboost_program_options
OUTPUT = solver.out

# Source files
SOURCES := $(SOURCEDIR)SolveMain.cpp

# Main target
all: $(OUTPUT)

# Create solver executable
$(OUTPUT): $(SOURCES)
	$(CC) $(FLAGS) $^ $(LIBS) -o $@

# Rule to clean up executable
clean:
	rm -f $(OUTPUT)

# Run the solver with mpiexec
run: $(OUTPUT)
	mpiexec -np 1 ./$(OUTPUT)