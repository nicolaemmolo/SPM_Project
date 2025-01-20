//
// Sequential and Parallel code of the first SPM Assignment
//		with Wavefront computation
//      MPI version + OMP
//
// compile:
// mpicxx -std=c++20 -O3 -fopenmp -I -Iinclude UTWavefrontMPIOMP.cpp -o UTWMPIOMP
//

#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <random>
#include <omp.h>


#ifndef PRINT_MESSAGE
    #define PRINT_MESSAGE 0
#endif

#ifndef PRINT_MATRIX
    #define PRINT_MATRIX 0
#endif

#ifndef PRINT_LAST_ELEMENT
	#define PRINT_LAST_ELEMENT 1
#endif


#define DEFAULT_DIM 3       // Default size of the matrix (NxN)
#define DEFAULT_NTHREADS 2	// Default number of threads
#define DEFAULT_NODES 2     // Default number of processes
#define DEFAULT_LOG_FILE "wavefront_results_MPIOMP.csv" // Default log file name

// Macro to calculate the index of the trinagular matrix element (element, diagonal, size)
#define INDEX(i,k,N) ((k)==0 ? (i) : ((N)*k - ((k)*(k-1))/2 + (i)))


// ---------------------- General functions ---------------------- //


/* Calculate and update a matrix element using the dot product
 * @param M: matrix
 * @param N: size of the matrix
 * @param i: element of the diagonal
 * @param k: diagonal index
 */
void compute_diagonal_element(std::vector<double> &M, const uint64_t &N, const uint64_t &i, const uint64_t &k) {
    double result = 0.0;

	// Calculate the dot product
    for (uint64_t j = 0; j < k; ++j) {
        result += M[INDEX(i, j, N)] * M[INDEX(k+i-j, j, N)];
    }

    M[INDEX(i, k, N)] = std::cbrt(result); // Update the element i for the diagonal k
}

/* Print matrix
 * @param M: matrix
 * @param N: size of the matrix
 */
void print_matrix(const std::vector<double> &M, uint64_t N) {
    std::cout << std::fixed << std::setprecision(2);

    for(uint64_t i = 0; i < N; ++i) {
		for (uint64_t j = 0; j < i; ++j) {
			std::cout << 0.00 << " ";  // Lower triangular part
		}
        for(uint64_t k = 0; k < N-i; ++k) {
            std::cout << M[INDEX(i,k,N)] << " ";  // Upper triangular part
		}
        std::cout << std::endl;
    }
}

/* Print M as a 1D array
 * @param M: matrix
 * @param total_elements: total number of elements in the upper triangular matrix
 */
void print_M(const std::vector<double> &M, uint64_t total_elements) {
	for (uint64_t i = 0; i < total_elements; ++i) {
		std::cout << M[i] << " ";
	}
	std::cout << std::endl;
}

/*
 * Print last element of the matrix (last diagonal element)
 * @param M: matrix
 * @param total_elements: total number of elements in the upper triangular matrix
 */
void print_last_element(const std::vector<double> &M, uint64_t total_elements) {
	std::cout << M[total_elements-1] << std::endl;
}


// ---------------------- Wavefront ---------------------- //


/* Wavefront (parallel MPI version)
 * @param M: matrix
 * @param N: size of the matrix
 */
void wavefront_parallel_mpi_omp(std::vector<double> &M, const uint64_t &N, const uint64_t &T, int rank, const int size) {
    for(int k = 1; k < N; ++k) {    // for each upper diagonal
        std::vector<int> counts(size);
        std::vector<int> displs(size);
     
        // Compute chunk_size and remainder given the k-th diagonal
        int chunk_size = (N-k) / size;
        int remainder = (N-k) % size;
        // Recompute the interval for each worker
        for (int i = 0; i < size; ++i){
            auto start = i * chunk_size + (i < remainder ? i : remainder);
            auto end = (i + 1) * chunk_size + (i < remainder ? (i + 1) : remainder);
            // Compute offsets and displacements
            counts[i] = end - start;
            displs[i] = start;
        }

        std::vector<double> buffer; // Buffer to hold the computed results
        std::vector<double> collect((N-k)); // Vector to collect all results from all processes

        // Process elements in the k-th diagonal, subdivided between workers
        #pragma omp parallel
        {
            std::vector<double> local_buffer; // Local buffer to hold the computed results
            #pragma omp for schedule(dynamic)
            for (int i = displs[rank]; i < (displs[rank] + counts[rank]); ++i) {
                compute_diagonal_element(M, N, i, k);
                local_buffer.push_back(M[INDEX(i, k, N)]);
            }

            // Synchronize threads and combine results
            #pragma omp critical
            buffer.insert(buffer.end(), local_buffer.begin(), local_buffer.end());
        }

        // Gather results from all processes
        MPI_Allgatherv(buffer.data(), counts[rank], MPI_DOUBLE, collect.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
   
        // Update diagonal elements
        for(int i = 0; i < (N-k); ++i){
            if (i >= displs[rank] && i < (displs[rank] + counts[rank]))
                continue;
            M[INDEX(i, k, N)] = collect[i];
        }
    }
}


/* ---------------------- Main function ---------------------- */


/* Main function
 * @param argc: number of arguments
 * @param argv: arguments
 * @return: 0 if successful
 */
int main(int argc, char *argv[]) {
    uint64_t N                = DEFAULT_DIM;
    uint64_t T                = DEFAULT_NODES;
    std::string log_file_name = DEFAULT_LOG_FILE;
    uint64_t nodes            = DEFAULT_NODES;

    // Initialize MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (provided < MPI_THREAD_FUNNELED) {
        if (rank == 0) {
            std::cerr << "MPI implementation does not provide full thread support!" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
        return -1;
    }

    
    // Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 3 && argc != 4) {
        if (rank == 0) {
            std::printf("use: %s [N] [log_file_name] [nodes]\n", argv[0]);
            std::printf("     N    : size of the square matrix\n");
            std::printf("     log_file_name : name of the log file\n");
            std::printf("     nodes: number of nodes\n");
            std::printf("     T    : number of threads\n");
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
                if (argc > 4) {
                    T = std::stol(argv[4]);
                }
            }
        }
    }

    omp_set_num_threads(T);

    // Allocate the matrix as a 1D array
    uint64_t total_elements = (N*(N+1))/2; // Total number of elements in the upper triangular matrix
    std::vector<double> M(total_elements, 0.0);

    if (rank == 0) {
        // Init function (initialize the diagonal elements)
        auto init = [&]() {
            for (uint64_t i = 0; i < N; ++i) {
                M[i] = static_cast<double>(i+1) / static_cast<double>(N);
            }
        };

        init();
    }

    // Broadcast the initialized matrix to all processes
    #pragma omp master // Only the master thread should execute this block
    MPI_Bcast(M.data(), total_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double execution_time = -1;

    // Parallel MPI
    if (rank == 0 && PRINT_MESSAGE) std::printf("------ Parallel MPI Execution ------\n");
    MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes before timing
    double start_time = MPI_Wtime();
    wavefront_parallel_mpi_omp(M, N, T, rank, size);
    double end_time = MPI_Wtime();
    execution_time = end_time - start_time;

    if (rank == 0) {
        if (PRINT_MATRIX) print_matrix(M, N);
        if (PRINT_MATRIX) print_M(M,total_elements);
	    if (PRINT_LAST_ELEMENT) print_last_element(M,total_elements);

        // Write the execution times to a file
        std::ofstream file;
        file.open(log_file_name, std::ios_base::app);
        file << N << "," << nodes << T << "," << execution_time << "\n";
        file.close();
        std::printf("Execution time: %f\n", execution_time);
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}