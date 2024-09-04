//
// Sequential and Parallel code of the first SPM Assignment
//		with Wavefront computation
//
// compile:
// mpicxx -std=c++20 -O3 -I -Iinclude UTWavefrontMPI_OMP.cpp -o UTWMPIOMP
//

#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <omp.h>
#include <mpi.h>


#define DEFAULT_DIM 3 		// Default size of the matrix (NxN)
#define DEFAULT_NTHREADS 2	// Default number of threads
#define DEFAULT_NODES 1     // Default number of nodes
#define DEFAULT_MODE "ps" 	// Default execution mode
#define DEFAULT_LOG_FILE "wavefront_results.csv"	// Default log file name


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


/* Wavefront (parallel version with static scheduling using OpenMP and MPI)
 * @param M: matrix
 * @param N: size of the matrix
 * @param T: number of threads
 * @param rank: MPI rank
 * @param size: number of MPI processes
 */
void wavefront_parallel_static_omp_mpi(std::vector<std::vector<double>> &M, const uint64_t &N, const uint32_t &T, const int &rank, const int &size) {
    #pragma omp parallel for num_threads(T) schedule(static)
    for (uint64_t k=1; k<N; ++k) { // For each upper diagonal
        for (uint64_t m=0; m<(N-k); ++m) { // For each element in the diagonal
            if (m % size == rank) { // Assign work based on rank
                compute_diagonal_element(M, N, m, k);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes
    }
}

/* Wavefront (parallel version with dynamic scheduling using OpenMP and MPI)
 * @param M: matrix
 * @param N: size of the matrix
 * @param T: number of threads
 * @param rank: MPI rank
 * @param size: number of MPI processes
 */
void wavefront_parallel_dynamic_omp_mpi(std::vector<std::vector<double>> &M, const uint64_t &N, const uint32_t &T, const int &rank, const int &size) {
    #pragma omp parallel for num_threads(T) schedule(dynamic)
    for (uint64_t k=1; k<N; ++k) { // For each upper diagonal
        for (uint64_t m=0; m<(N-k); ++m) { // For each element in the diagonal
            if (m % size == rank) { // Assign work based on rank
                compute_diagonal_element(M, N, m, k);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes
    }
}

int main(int argc, char *argv[]) {
    uint64_t N                = DEFAULT_DIM;
    uint32_t T                = DEFAULT_NTHREADS;
    std::string mode          = "ps"; // Default mode
    std::string log_file_name = DEFAULT_LOG_FILE;
    uint64_t nodes            = DEFAULT_NODES;

    # Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 3 && argc != 4 && argc != 5) {
        if (rank == 0) {
            std::printf("use: %s [N] [T] [mode] [log_file_name] [nodes]\n", argv[0]);
            std::printf("     N    : size of the square matrix\n");
            std::printf("     T    : number of threads\n");
            std::printf("     mode : execution mode ('ps' for parallel static, 'pd' for parallel dynamic)\n");
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
            T = std::stol(argv[2]);
            if (argc > 3) {
                mode = argv[3];
                if (argc > 4) {
                    log_file_name = argv[4];
                    if (argc > 5) {
                        nodes = std::stol(argv[5]);
                    }
                }
            }
        }
    }

    // Allocate the matrix
    std::vector<std::vector<double>> M(N, std::vector<double>(N, 0.0));
    
    // Init function
    auto init = [&]() {
        for (uint64_t m = 0; m < N; ++m) {
            M[m][m] = static_cast<double>(m + 1) / static_cast<double>(N);
        }
    };
    
    init();

    double execution_time = -1;

    // Parallel static
    if (mode == "ps") {
        if (rank == 0 && PRINT_MESSAGE) std::printf("------ Parallel Static Execution ------\n");
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes before timing
        double start_time = MPI_Wtime();
        wavefront_parallel_static_omp_mpi(M, N, T, rank, size);
        double end_time = MPI_Wtime();
        execution_time = end_time - start_time;
    }

    // Parallel dynamic
    if (mode == "pd") {
        if (rank == 0 && PRINT_MESSAGE) std::printf("------ Parallel Dynamic Execution ------\n");
        MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes before timing
        double start_time = MPI_Wtime();
        wavefront_parallel_dynamic_omp_mpi(M, N, T, rank, size);
        double end_time = MPI_Wtime();
        execution_time = end_time - start_time;
    }


    // Write the execution times to a file
    if (rank == 0) {
        if (PRINT_MATRIX) print_matrix(M, N);

        std::ofstream file;
        file.open(log_file_name, std::ios_base::app);
        file << N << "," << T << "," << nodes << "," << mode << "," << execution_time << "\n";
        file.close();
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
