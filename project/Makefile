# Compilation variables
CXX = g++
CXXFLAGS = -std=c++20 -O3 -march=native -Iinclude

MPI_CXX = mpicxx
MPI_CXXFLAGS = -std=c++20 -O3 -I -Iinclude

# Executable names
EXEC_UTW = UTW
EXEC_UTWFF = UTWFF
EXEC_UTWOMP = UTWOMP
EXEC_UTWMPI = UTWMPI

# Source files
SRC_UTW = UTWavefront.cpp
SRC_UTWFF = UTWavefrontFastFlow.cpp
SRC_UTWOMP = UTWavefrontOpenMP.cpp
SRC_UTWMPI = UTWavefrontMPI.cpp

# Default rule: compile all executables
all: $(EXEC_UTW) $(EXEC_UTWFF) $(EXEC_UTWOMP) $(EXEC_UTWMPI)

# Compile the UTW executable
$(EXEC_UTW): $(SRC_UTW)
	$(CXX) $(CXXFLAGS) $< -o $@

# Compile the UTWFF executable
$(EXEC_UTWFF): $(SRC_UTWFF)
	$(CXX) $(CXXFLAGS) -Iinclude/fastflow-master $< -o $@

# Compile the UTWOMP executable
$(EXEC_UTWOMP): $(SRC_UTWOMP)
	$(CXX) $(CXXFLAGS) -fopenmp $< -o $@

# Compile the UTWMPI executable
$(EXEC_UTWMPI): $(SRC_UTWMPI)
	$(MPI_CXX) $(MPI_CXXFLAGS) $< -o $@

# Rule to clean up executables and object files
clean:
	rm -f $(EXEC_UTW) $(EXEC_UTWFF) $(EXEC_UTWOMP) $(EXEC_UTWMPI)

# Add other rules if necessary
# example: test, install, etc.

.PHONY: all clean
