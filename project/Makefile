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
EXEC_UTWMPIOMP = UTWMPIOMP

# Source files
SRC_UTW = UTWavefront.cpp
SRC_UTWFF = UTWavefrontFastFlow.cpp
SRC_UTWOMP = UTWavefrontOpenMP.cpp
SRC_UTWMPI = UTWavefrontMPI.cpp
SRC_UTWMPIOMP = UTWavefrontMPI_OMP.cpp

# Default rule: compile all executables
all: $(EXEC_UTW) $(EXEC_UTWFF) $(EXEC_UTWOMP) $(EXEC_UTWMPI) $(EXEC_UTWMPIOMP)

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

# Compile the UTWMPIOMP executable
$(EXEC_UTWMPIOMP): $(SRC_UTWMPIOMP)
	$(MPI_CXX) $(MPI_CXXFLAGS) -fopenmp $< -o $@

# Rule to clean up executables and object files
clean:
	rm -f $(EXEC_UTW) $(EXEC_UTWFF) $(EXEC_UTWOMP) $(EXEC_UTWMPI) $(EXEC_UTWMPIOMP)

# Declare the rules as phony
.PHONY: all clean
