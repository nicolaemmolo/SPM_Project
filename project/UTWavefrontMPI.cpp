//
// Sequential and Parallel code of the first SPM Assignment
//		with Wavefront computation
//      MPI version
//
// compile:
// mpicxx -std=c++20 -O3 -I -Iinclude UTWavefrontMPI.cpp -o UTWMPI
//

#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <mpi.h>


#ifndef PRINT_MESSAGE
    #define PRINT_MESSAGE 0
#endif

#ifndef PRINT_MATRIX
    #define PRINT_MATRIX 1
#endif

#define DEFAULT_DIM 3       // Default size of the matrix (NxN)
#define DEFAULT_NODES 2     // Default number of threads
#define DEFAULT_LOG_FILE "wavefront_results_MPI.csv" // Default log file name


/* Calculate dot product of two vectors
 * @param v1: first vector
 * @param v2: second vector
 * @return: dot product of the two vectors
 */
double dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

/* Calculate and update a matrix element using the dot product of two vectors
 * @param M: matrix
 * @param N: size of the matrix
 * @param m: row index
 * @param k: column index
 */
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
}

/* Print matrix
 * @param M: matrix
 * @param N: size of the matrix
 */
void print_matrix(const std::vector<std::vector<double>> &M, uint64_t N) {
    std::cout << std::fixed << std::setprecision(2);
    for (uint64_t i=0; i<N; ++i) {
        for (uint64_t j=0; j<N; ++j) {
            std::cout << M[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


/* Wavefront (parallel MPI version)
 * @param M: matrix
 * @param N: size of the matrix
 */
void wavefront_parallel_mpi(std::vector<std::vector<double>> &M, const uint64_t &N) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    for (uint64_t k=1; k<N; ++k) { // For each upper diagonal
        for (uint64_t m=0; m<(N-k); ++m) { // For each element in the diagonal
            if (m % size == rank) { // Assign work based on rank
                compute_diagonal_element(M, N, m, k);
            }
            MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes
        }
    }
}


/* Main function
 * @param argc: number of arguments
 * @param argv: arguments
 * @return: 0 if successful
 */
int main(int argc, char *argv[]) {
    uint64_t N                = DEFAULT_DIM;
    std::string log_file_name = DEFAULT_LOG_FILE;
    uint64_t nodes            = DEFAULT_NODES;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 3 && argc != 4) {
        if (rank == 0) {
            std::printf("use: %s [N] [log_file_name] [nodes]\n", argv[0]);
            std::printf("     N    : size of the square matrix\n");
            std::printf("     log_file_name : name of the log file\n");
            std::printf("     nodes: number of nodes\n");
        }
        MPI_Finalize();
        return -1;
    }

    // Set the parameters (if any)
    if (argc > 1) {
        N = std::stol(argv[1]);
        if (argc > 2) {
            log_file_name = argv[2];
            if (argc > 3) {
                nodes = std::stol(argv[3]);
            }
        }
    }

    // Allocate the matrix
    std::vector<std::vector<double>> M(N, std::vector<double>(N, 0.0));

    // Init function
    auto init = [&]() {
        for (uint64_t m=0; m<N; ++m) {
            M[m][m] = static_cast<double>(m+1) / static_cast<double>(N);
        }
    };
    
    init();

    double execution_time = -1;

    // Parallel MPI
    if (rank == 0 && PRINT_MESSAGE) std::printf("------ Parallel MPI Execution ------\n");
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes before timing
    double start_time = MPI_Wtime();
    wavefront_parallel_mpi(M, N);
    double end_time = MPI_Wtime();
    execution_time = end_time - start_time;

    if (rank == 0) {
        if (PRINT_MATRIX) print_matrix(M, N);

        // Write the execution times to a file
        std::ofstream file;
        file.open(log_file_name, std::ios_base::app);
        file << N << "," << nodes <<"," << execution_time << "\n";
        file.close();
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}