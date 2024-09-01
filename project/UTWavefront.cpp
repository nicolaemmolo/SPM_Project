//
// Sequential and Parallel code of the first SPM Assignment
//		with Wavefront computation
//		C++ Thread version
//
// compile:
// g++ -std=c++20 -O3 -march=native -Iinclude UTWavefront.cpp -o UTW
//

#include <iostream>
#include <vector>
#include <thread>
#include <random>
#include <cassert>
#include <barrier>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <hpc_helpers.hpp>
#include <threadPool.hpp>

#ifndef PRINT_MESSAGE
	#define PRINT_MESSAGE 0
#endif

#ifndef PRINT_MATRIX
	#define PRINT_MATRIX 0
#endif

#define DEFAULT_DIM 3 		// Default size of the matrix (NxN)
#define DEFAULT_NTHREADS 2	// Default number of threads
#define DEFAULT_MODE "s" 	// Default execution mode
#define DEFAULT_LOG_FILE "wavefront_results_UTW.csv"	// Default log file name


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


/* Wavefront (sequential version)
 * @param M: matrix
 * @param N: size of the matrix
 */
void wavefront_sequential(std::vector<std::vector<double>> &M, const uint64_t &N) {
    for (uint64_t k=1; k<N; ++k) { // For each upper diagonal
        for (uint64_t m=0; m<(N-k); ++m) { // For each element in the diagonal
            compute_diagonal_element(M, N, m, k);
        }
    }
}


/* Wavefront (parallel version with static scheduling)
 * @param M: matrix
 * @param N: size of the matrix
 * @param T: number of threads
 */
void wavefront_parallel_static(std::vector<std::vector<double>> &M, const uint64_t &N, const uint64_t &T) {
	std::barrier barrier(T);

	auto static_task = [&] (const uint64_t id) -> void {
		for (uint64_t k=1; k<N; ++k) { // For each upper diagonal
			if (id >= N-k) { // If the thread is not needed
				barrier.arrive_and_drop();
				return;
			}
			for (uint64_t i=id; i<N-k; i+=T) { // For each assigned element in the diagonal 
				compute_diagonal_element(M, N, i, k);
			}
			barrier.arrive_and_wait();
		}
	};

	std::vector<std::thread> threads;
	for (uint64_t id=0; id<T; ++id)
		threads.emplace_back(static_task, id);

	for (auto &thread : threads)
		thread.join();
}


/* Wavefront (parallel version with dynamic scheduling)
 * @param M: matrix
 * @param N: size of the matrix
 * @param T: number of threads
 */
void wavefront_parallel_dynamic(std::vector<std::vector<double>> &M, const uint64_t &N, const uint64_t &T) {
    /*
	ThreadPool pool(T);

    for (uint64_t k=1; k<N; ++k) { // For each upper diagonal
        uint64_t active_threads = std::min(T, N-k); // Calculate the active threads for this diagonal
        std::barrier barrier(active_threads); // Create a new barrier for this iteration

        if ((N-k) < T) { // If the diagonal is smaller than the number of threads
            for (uint64_t i=(N-k); i<T; ++i) { // For each extra thread
                pool.enqueue([&barrier]() {
                    barrier.arrive_and_wait();
                });
            }
        }

        for (uint64_t i=0; i<(N-k); ++i) { // For each element in the diagonal
            bool block = (i >= (N-k-T) || (N-k) < T);
            pool.enqueue([&barrier, i, k, &M, N, block]() { 
                compute_diagonal_element(M, N, i, k);
                if (block) { // If the thread should wait
                    barrier.arrive_and_wait();
                }
            });
        }
    } */


	std::barrier bar(T);

	auto process_element = [&](uint64_t index, bool block) {
		uint64_t i = index / N;
		uint64_t k = index % N - i;
		compute_diagonal_element(M, N, i, k);
		if (block) { // If the thread should wait
			bar.arrive_and_wait();
		}
	};

	auto wait = [&] () {
		bar.arrive_and_wait();
	};

	ThreadPool TP(T);
	for (uint64_t k=0; k<N; ++k) { // For each upper diagonal
		if ((N-k)<T) { // If the diagonal is smaller than the number of threads
			for (uint64_t i=(N-k); i<T; ++i) { // For each extra thread
				TP.enqueue(wait);
			}
		}
		for (uint64_t i=0; i<(N-k); ++i) { // For each element in the diagonal
			bool block = (i>=(N-k-T) || (N-k)<T) ? true : false;
			TP.enqueue(process_element, i*N+(i+k), block);
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
	uint64_t T                = DEFAULT_NTHREADS;
	std::string mode          = DEFAULT_MODE;
	std::string log_file_name = DEFAULT_LOG_FILE;

	// Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 3 && argc != 4 && argc != 5) {
        std::printf("use: %s [N] [T] [mode] [log_file_name]\n", argv[0]);
        std::printf("     N    : size of the square matrix\n");
        std::printf("     T    : number of threads\n");
        std::printf("     mode : execution mode ('s' for sequential, 'ps' for parallel static, 'pd' for parallel dynamic)\n");
		std::printf("     log_file_name : name of the log file\n");
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
				}
			}
		}
	}

	// Allocate the matrix
	std::vector<std::vector<double>> M(N, std::vector<double>(N, 0.0));

	// Init function
	auto init=[&]() {
		for (uint64_t m=0; m<N; ++m) {
			M[m][m] = static_cast<double>(m+1) / static_cast<double>(N);
		}
	};
	
	init();

	double execution_time=-1;

	// Sequential
	if (mode == "s") {
		if (PRINT_MESSAGE) std::printf("------ Sequential Execution ------\n");
		TIMERSTART(wavefront_sequential);
		wavefront_sequential(M, N); 
		TIMERSTOP(wavefront_sequential, execution_time);
	}
	
	// Parallel static
	if (mode == "ps") {
		if (PRINT_MESSAGE) std::printf("------ Parallel Static Execution ------\n");
		TIMERSTART(wavefront_parallel_static);
		wavefront_parallel_static(M, N, T); 
		TIMERSTOP(wavefront_parallel_static, execution_time);
	}

	// Parallel dynamic
	if (mode == "pd") {
		if (PRINT_MESSAGE) std::printf("------ Parallel Dynamic Execution ------\n");
		TIMERSTART(wavefront_parallel_dynamic);
		wavefront_parallel_dynamic(M, N, T); 
		TIMERSTOP(wavefront_parallel_dynamic, execution_time);
	}

	if (PRINT_MATRIX) print_matrix(M,N);

	
	// Write the execution times to a file
	std::ofstream file;
	file.open(log_file_name, std::ios_base::app);
	file << N << "," << T << "," << mode << "," << execution_time << "\n";
	file.close();

	return 0;
}