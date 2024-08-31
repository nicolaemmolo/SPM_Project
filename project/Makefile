# Variabili di compilazione
CXX = g++
CXXFLAGS = -std=c++20 -O3 -march=native -Iinclude

# Nomi degli eseguibili
EXEC_UTW = UTW
EXEC_UTWFF = UTWFF
EXEC_UTWOMP = UTWOMP

# File sorgente
SRC_UTW = UTWavefront.cpp
SRC_UTWFF = UTWavefrontFastFlow.cpp
SRC_UTWOMP = UTWavefrontOpenMP.cpp

# Regola di default: compila tutti gli eseguibili
all: $(EXEC_UTW) $(EXEC_UTWFF) $(EXEC_UTWOMP)

# Compilazione dell'eseguibile UTW
$(EXEC_UTW): $(SRC_UTW)
	$(CXX) $(CXXFLAGS) $< -o $@

# Compilazione dell'eseguibile UTWFF
$(EXEC_UTWFF): $(SRC_UTWFF)
	$(CXX) $(CXXFLAGS) -Iinclude/fastflow-master $< -o $@

# Compilazione dell'eseguibile UTWOMP
$(EXEC_UTWOMP): $(SRC_UTWOMP)
	$(CXX) $(CXXFLAGS) -fopenmp $< -o $@

# Regola per pulire gli eseguibili e i file oggetto
clean:
	rm -f $(EXEC_UTW) $(EXEC_UTWFF) $(EXEC_UTWOMP)

# Aggiungi altre regole se necessario
# esempio: test, install, ecc.

.PHONY: all clean
