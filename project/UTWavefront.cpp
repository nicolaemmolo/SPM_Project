//
// Sequential and Parallel code of the first SPM Assignment
//		with Wavefront computation
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
	#define PRINT_MATRIX 1
#endif

#define DEFAULT_DIM 3 		// default size of the matrix (NxN)
#define DEFAULT_NTHREADS 2	// default number of threads
#define DEFAULT_LOG_FILE "wavefront_results.csv"	// default log file name


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


// wavefront (sequential version)
void wavefront_sequential(std::vector<std::vector<double>> &M, const uint64_t &N) {
    for (uint64_t k=1; k<N; ++k) { // for each upper diagonal
        for (uint64_t m=0; m<(N-k); ++m) { // for each element in the diagonal
            compute_diagonal_element(M, N, m, k);
        }
    }
}


// wavefront (parallel version with static scheduling)
void wavefront_parallel_static(std::vector<std::vector<double>> &M, const uint64_t &N, const uint32_t &T) {

	std::barrier bar(T);

	auto static_task = [&] (const uint64_t id) -> void {
		for (uint64_t k=1; k<N; ++k) { // for each upper diagonal
			if (id >= N-k) {
				bar.arrive_and_drop();
				if (PRINT_MESSAGE) std::printf("Thread %lu: exit\n", id);
				return;
			}
			for (uint64_t i=id; i<N-k; i+=T) { // for each assigned element
				compute_diagonal_element(M, N, i, k);
				if (PRINT_MESSAGE) std::printf("Thread %lu: computed index %lu\n", id, i*N+(i+k));
			}
			bar.arrive_and_wait();
			if (PRINT_MESSAGE) std::printf("Thread %lu: unlocked from waiting\n", id);
		}
	};

	std::vector<std::thread> threads;
	for (uint64_t id=0; id<T; id++)
		threads.emplace_back(static_task, id);

	for (auto &thread : threads)
		thread.join();
}


// wavefront (parallel version with dynamic scheduling)
void wavefront_parallel_dynamic(std::vector<std::vector<double>> &M, const uint64_t &N, const uint32_t &T) {

	std::barrier bar(T);

	auto process_element = [&](uint64_t index, bool block) {
		uint64_t i = index / N;
		uint64_t k = index % N - i;
		compute_diagonal_element(M, N, i, k);
		if (block) {
			bar.arrive_and_wait();
			if (PRINT_MESSAGE) std::printf("Computed index %lu and waited\n", index);
		}
		else {
			if (PRINT_MESSAGE) std::printf("Computed index %lu without waiting\n", index);
		}
	};

	auto wait = [&] () {
		bar.arrive_and_wait();
		if (PRINT_MESSAGE) std::printf("Waited\n");
	};

	ThreadPool TP(T);
	for (uint64_t k=1; k<N; ++k) { // for each upper diagonal
		if ((N-k)<T) { // if the diagonal is smaller than the number of threads
			for (uint64_t i=(N-k); i<T; ++i) { // for each extra thread
				TP.enqueue(wait);
			}
		}
		for (uint64_t i=0; i<(N-k); ++i) { // for each elem. in the diagonal
			bool block = (i>=(N-k-T) || (N-k)<T) ? true : false;
			TP.enqueue(process_element, i*N+(i+k), block);
		}
	}
}


int main(int argc, char *argv[]) {
	uint64_t N                = DEFAULT_DIM;
	uint32_t T                = DEFAULT_NTHREADS;
	std::string log_file_name = DEFAULT_LOG_FILE;

	// Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 4 && argc != 5) {
        std::printf("use: %s [N] [T]\n", argv[0]);
        std::printf("     N size of the square matrix\n");
        std::printf("     T number of threads\n");
        return -1;
    }

	// Set the parameters (if any)
    if (argc > 1) {
        N = std::stol(argv[1]);
        if (argc > 2) {
            T = std::stol(argv[2]);
        }
    }

	// allocate the matrix
	std::vector<std::vector<double>> M(N, std::vector<double>(N, 0.0));


	uint64_t expected_totaltime=0;
	// init function
	auto init=[&]() {
		for (uint64_t m=0; m<N; ++m) {
			M[m][m] = static_cast<double>(m+1) / static_cast<double>(N);
		}
	};
	
	init();

	// sequential
	double sequential_time=-1;
	if (PRINT_MESSAGE) std::printf("------ Sequential Execution ------\n");
	TIMERSTART(wavefront_sequential);
	wavefront_sequential(M, N); 
	TIMERSTOP(wavefront_sequential, sequential_time);
	
	if (PRINT_MATRIX) print_matrix(M,N);

	// parallel static
	double parallel_static_time=-1;
	if (PRINT_MESSAGE) std::printf("------ Parallel Static Execution ------\n");
	TIMERSTART(wavefront_parallel_static);
	wavefront_parallel_static(M, N, T); 
	TIMERSTOP(wavefront_parallel_static, parallel_static_time);

	if (PRINT_MATRIX) print_matrix(M,N);

	// parallel dynamic
	double parallel_dynamic_time=-1;
	if (PRINT_MESSAGE) std::printf("------ Parallel Dynamic Execution ------\n");
	TIMERSTART(wavefront_parallel_dynamic);
	wavefront_parallel_dynamic(M, N, T); 
	TIMERSTOP(wavefront_parallel_dynamic, parallel_dynamic_time);

	if (PRINT_MATRIX) print_matrix(M,N);

	// write the execution times to a file
	std::ofstream file;
	file.open(log_file_name, std::ios_base::app);
	file << N << "," << T << ","
		<< expected_totaltime/1000000 << "," << sequential_time << ","
		<< parallel_dynamic_time << "," << parallel_static_time << "\n";
	file.close();

	return 0;
}