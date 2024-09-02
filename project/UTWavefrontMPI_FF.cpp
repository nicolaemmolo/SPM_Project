//
// Sequential and Parallel code of the first SPM Assignment
//		with Wavefront computation
//
// compile:
// g++ -std=c++20 -O3 -march=native -Iinclude UTWavefront.cpp -o UTW
//

#include <mpi.h>
#include <iostream>
#include <vector>
#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <hpc_helpers.hpp>

using namespace ff;

#ifndef PRINT_MESSAGE
    #define PRINT_MESSAGE 0
#endif

#ifndef PRINT_MATRIX
    #define PRINT_MATRIX 1
#endif

#define DEFAULT_DIM 3         // default size of the matrix (NxN)
#define DEFAULT_NTHREADS 2    // default number of threads
#define DEFAULT_MODE "pd"     // default execution mode
#define DEFAULT_LOG_FILE "wavefront_results.csv"    // default log file name

// Calculate dot product of two vectors
double dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

// Calculate and update matrix element with dot product of two vectors
void compute_diagonal_element(std::vector<std::vector<double>> &M, const uint64_t &N, const uint64_t &m, const uint64_t &k) {
    std::vector<double> row_vector(k);
    std::vector<double> col_vector(k);
    
    // Fill vectors with the corresponding elements from row m and column m+k
    for (uint64_t i=0; i<k; ++i) {
        row_vector[i] = M[m][m+i];
        col_vector[i] = M[m+i+1][m+k];
    }
    
    // Compute the dot product
    M[m][m+k] = dot_product(row_vector, col_vector);
    std::printf("compute\n");
}

// Print matrix
void print_matrix(const std::vector<std::vector<double>> &M, uint64_t N) {
    std::cout << std::fixed << std::setprecision(2);
    for (uint64_t i=0; i<N; ++i) {
        for (uint64_t j=0; j<N; ++j) {
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// wavefront (parallel version with dynamic scheduling using FastFlow)
void wavefront_parallel_dynamic_ff(std::vector<std::vector<double>> &M, const uint64_t &N, const uint32_t &T) {
    ParallelFor pf(T);
    
    for (uint64_t k=1; k<N; ++k) { // for each upper diagonal
        pf.parallel_for(0, N-k, [&](const long i) {
            compute_diagonal_element(M, N, i, k);
        });
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint64_t N = DEFAULT_DIM;
    uint32_t T = DEFAULT_NTHREADS;
    std::string mode = DEFAULT_MODE;
    std::string log_file_name = DEFAULT_LOG_FILE;

    // Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 3 && argc != 4) {
        if (rank == 0) {
            std::printf("use: %s [N] [T] [mode]\n", argv[0]);
            std::printf("     N    : size of the square matrix\n");
            std::printf("     T    : number of threads\n");
            std::printf("     mode : execution mode ('s' for sequential, 'ps' for parallel static, 'pd' for parallel dynamic)\n");
        }
        MPI_Finalize();
        return -1;
    }

    // Set the parameters (if any)
    if (argc > 1) {
        N = std::stol(argv[1]);
        if (argc > 2) {
            T = std::stol(argv[2]);
            if (argc > 3) {
                mode = argv[3];
            }
        }
    }

    // Determine the submatrix each process will handle
    uint64_t rows_per_proc = N / size;
    uint64_t start_row = rank * rows_per_proc;
    uint64_t end_row = (rank == size - 1) ? N : start_row + rows_per_proc;

    // Allocate the local submatrix
    std::vector<std::vector<double>> M_local(rows_per_proc, std::vector<double>(N, 0.0));

    // Initialize the local submatrix
    for (uint64_t m = 0; m < rows_per_proc; ++m) {
        M_local[m][m + start_row] = static_cast<double>(m + 1 + start_row) / static_cast<double>(N);
    }

    double execution_time = -1;

    // Perform the wavefront computation
    if (mode == "pd") {
        if (PRINT_MESSAGE && rank == 0) std::printf("------ Parallel Dynamic Execution with MPI ------\n");
        TIMERSTART(wavefront_parallel_dynamic_ff);
        wavefront_parallel_dynamic_ff(M_local, N, T);
        TIMERSTOP(wavefront_parallel_dynamic_ff, execution_time);
    }

    // Gather the results in the root process
    std::vector<std::vector<double>> M;
    if (rank == 0) {
        M.resize(N, std::vector<double>(N, 0.0));
    }

    for (uint64_t i = 0; i < rows_per_proc; ++i) {
        MPI_Gather(M_local[i].data(), N, MPI_DOUBLE, (rank == 0) ? M[i + start_row].data() : nullptr, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    if (rank == 0 && PRINT_MATRIX) print_matrix(M, N);

    // Write the execution times to a file
    if (rank == 0) {
        std::ofstream file;
        file.open(log_file_name, std::ios_base::app);
        file << N << "," << T << "," << mode << "," << execution_time << "\n";
        file.close();
    }

    MPI_Finalize();
    return 0;
}
