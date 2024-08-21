//
// Sequential and Parallel code of the first SPM Assignment
//		with Wavefront computation
//
// compile:
// g++ -std=c++20 -O3 -march=native -Iinclude UTWavefrontMPI.cpp -o UTWMPI
//

#include <iostream>
#include <vector>
#include <mpi.h>
#include <random>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <cassert>

#define DEFAULT_DIM 3 		// default size of the matrix (NxN)
#define DEFAULT_LOG_FILE "wavefront_results_mpi.csv"	// default log file name

// Funzione di prodotto scalare
double dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

// Funzione per calcolare e aggiornare un elemento della matrice
void compute_diagonal_element(std::vector<std::vector<double>> &M, const uint64_t &N, const uint64_t &m, const uint64_t &k) {
    std::vector<double> row_vector(k);
    std::vector<double> col_vector(k);
    
    for (uint64_t i = 0; i < k; ++i) {
        row_vector[i] = M[m][m + i];
        col_vector[i] = M[m + i + 1][m + k];
    }
    
    M[m][m + k] = dot_product(row_vector, col_vector);
}

// Funzione per stampare la matrice
void print_matrix(const std::vector<std::vector<double>> &M, uint64_t N) {
    std::cout << std::fixed << std::setprecision(2);
    for (uint64_t i = 0; i < N; ++i) {
        for (uint64_t j = 0; j < N; ++j) {
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint64_t N = DEFAULT_DIM;
    std::string log_file_name = DEFAULT_LOG_FILE;

    if (argc > 1) {
        N = std::stol(argv[1]);
    }

    std::vector<std::vector<double>> M(N, std::vector<double>(N, 0.0));

    // Funzione di inizializzazione della matrice
    auto init = [&]() {
        for (uint64_t m = 0; m < N; ++m) {
            M[m][m] = static_cast<double>(m + 1) / static_cast<double>(N);
        }
    };

    init();

    double execution_time = -1.0;
    double start_time = MPI_Wtime();

    for (uint64_t k = 1; k < N; ++k) { // Per ogni diagonale superiore
        for (uint64_t m = rank; m < (N - k); m += size) { // Ogni processo si occupa di parte della diagonale
            compute_diagonal_element(M, N, m, k);
            if (m + 1 < (N - k)) {
                MPI_Send(&M[m][m + k], 1, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD);
            }
        }

        if (rank > 0) {
            MPI_Recv(&M[(rank - 1) % size][(rank - 1) % size + k], 1, MPI_DOUBLE, (rank - 1) % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    double end_time = MPI_Wtime();
    execution_time = end_time - start_time;

    if (rank == 0) {
        if (argc > 2 && std::string(argv[2]) == "print") {
            print_matrix(M, N);
        }

        std::ofstream file;
        file.open(log_file_name, std::ios_base::app);
        file << N << "," << size << "," << execution_time << "\n";
        file.close();
    }

    MPI_Finalize();
    return 0;
}
