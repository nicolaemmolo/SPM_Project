//
// Sequential and Parallel code of the first SPM Assignment
// 		with Wavefront computation
//		FastFlow version
//
// compile:
// g++ -std=c++20 -O3 -march=native -Iinclude -Iinclude/fastflow-master UTWavefrontFastFlow.cpp -o UTWFF
//

#include <iostream>
#include <vector>
#include <random>
#include <cassert>
#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>
#include <ff/pipeline.hpp>
#include <ff/farm.hpp>
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
#define DEFAULT_MODE "pd" 	// default execution mode
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


// wavefront (parallel version with static scheduling using FastFlow)
void wavefront_parallel_static_ff(std::vector<std::vector<double>> &M, const uint64_t &N, const uint32_t &T) {
    ParallelFor pf(T);
    
	for (uint64_t k=1; k<N; ++k) { // for each upper diagonal
        pf.parallel_for(0, N-k, 1, 0, [&](const long i) {
            compute_diagonal_element(M, N, i, k);
        });
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

    /*
    ParallelFor pf(T);
    
    for (uint64_t k=1; k<N; ++k) { // for each upper diagonal
        pf.parallel_for(0, N-k, [&](const long i) {
            // Calculate the indices corresponding to the element
            uint64_t index = i * N + (i+k);
            uint64_t m = index / N;
            uint64_t n = index % N - m;
            
            compute_diagonal_element(M, N, m, n);

            if (i >= (N-k-T) || (N-k) < T) {
                pf.parallel_for(0, 1, [&](const long) {
                    // Barrier synchronization point (dummy operation for waiting)
                    std::vector<double> dummy(1000000, 0.0);
                });
            }
        });
    }*/
}

int main(int argc, char *argv[]) {
    uint64_t N                = DEFAULT_DIM;
    uint32_t T                = DEFAULT_NTHREADS;
    std::string mode          = DEFAULT_MODE;
    std::string log_file_name = DEFAULT_LOG_FILE;

    // Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 3 && argc != 4) {
        std::printf("use: %s [N] [T] [mode]\n", argv[0]);
        std::printf("     N    : size of the square matrix\n");
        std::printf("     T    : number of threads\n");
        std::printf("     mode : execution mode ('s' for sequential, 'ps' for parallel static, 'pd' for parallel dynamic)\n");
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

    // allocate the matrix
    std::vector<std::vector<double>> M(N, std::vector<double>(N, 0.0));
    
    // init function
    auto init=[&]() {
        for (uint64_t m=0; m<N; ++m) {
            M[m][m] = static_cast<double>(m+1) / static_cast<double>(N);
        }
    };
    
    init();

    double execution_time=-1;

    // parallel static
	if (mode == "ps") {
		if (PRINT_MESSAGE) std::printf("------ Parallel Static Execution ------\n");
		TIMERSTART(wavefront_parallel_static_ff);
		wavefront_parallel_static_ff(M, N, T); 
		TIMERSTOP(wavefront_parallel_static_ff, execution_time);
	}

	// parallel dynamic
	if (mode == "pd") {
		if (PRINT_MESSAGE) std::printf("------ Parallel Dynamic Execution ------\n");
		TIMERSTART(wavefront_parallel_dynamic_ff);
		wavefront_parallel_dynamic_ff(M, N, T); 
		TIMERSTOP(wavefront_parallel_dynamic_ff, execution_time);
	}

	if (PRINT_MATRIX) print_matrix(M,N);
    

    // write the execution times to a file
    std::ofstream file;
    file.open(log_file_name, std::ios_base::app);
    file << N << "," << T << "," << mode << "," << execution_time << "\n";
    file.close();

    return 0;
}
