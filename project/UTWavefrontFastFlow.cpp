//
// Sequential code of the first SPM Assignment a.a. 23/24.
//
// compile:
// g++ -std=c++20 -O3 -march=native -Iinclude UTWavefrontFastFlow.cpp -o UTW
//
#include <iostream>
#include <vector>
#include <thread>
#include <random>
#include <cassert>
#include <barrier>
#include <fstream>
#include <numeric>
#include <hpc_helpers.hpp>
#include <threadPool.hpp>

#ifndef PRINT_MESSAGES
	#define PRINT_MESSAGES 1
#endif

#define DEFAULT_MIN 0 		// default minimum time (in microseconds)
#define DEFAULT_MAX 10		// default maximum time (in microseconds)
#define DEFAULT_DIM 10 		// default size of the matrix (NxN)
#define DEFAULT_NTHREADS 2	// default number of threads
#define DEFAULT_LOG_FILE "wavefront_results.csv"	// default log file name


int random(const int &min, const int &max) {
	static std::mt19937 generator(117);
	std::uniform_int_distribution<int> distribution(min,max);
	return distribution(generator);
};		

// Calculate dot product of two vectors
double dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

// Calculate and update matrix element with dot product of two vectors
void calculate_diagonal_element(const std::vector<double>& M, uint64_t i, uint64_t k, uint64_t N) {
	std::vector<double> v_m(M.begin() + i * N, M.begin() + i * N + k);
	std::vector<double> v_mk(M.begin() + (i + k) * N, M.begin() + (i + k) * N + k);
	M[i * N + (i+k)] = dot_product(v_m, v_mk);
}

// wavefront (sequential version)
void wavefront_sequential(const std::vector<double> &M, const uint64_t &N) {
	for(uint64_t k=0; k<N; ++k) { // for each upper diagonal
		for(uint64_t i=0; i<(N-k); ++i) { // for each element in the diagonal
			calculate_diagonal_element(M,i,k,N)
		}
	}
}

// wavefront (parallel version with static scheduling)
void wavefront_parallel_static(
	const std::vector<double> &M,
	const uint64_t &N,
	const uint32_t &T
	) {

		std::barrier bar(T);

		auto static_task = [&] (const uint64_t id) -> void {
			for (uint64_t k=0; k<N; ++k) { // for each upper diagonal
				if (id >= N-k) {
					bar.arrive_and_drop();
					if (PRINT_MESSAGES) std::printf("Thread %lu: exit\n", id);
					return;
				}
				for (uint64_t i=id; i<N-k; i+=T) { // for each assigned element
                	calculate_dot_product(M, i, k, N);
					if (PRINT_MESSAGES) std::printf("Thread %lu: computed index %lu\n", id, i*N+(i+k));
				}
				bar.arrive_and_wait();
				if (PRINT_MESSAGES) std::printf("Thread %lu: unlocked from waiting\n", id);
			}
		};

		std::vector<std::thread> threads;
		for (uint64_t id=0; id<T; id++)
			threads.emplace_back(static_task, id);

		for (auto &thread : threads)
			thread.join();
	}

// wavefront (parallel version with dynamic scheduling)
void wavefront_parallel_dynamic(
	const std::vector<double> &M,
	const uint64_t &N,
	const uint32_t &T
	) {

		std::barrier bar(T);

		auto process_element = [&](uint64_t index, bool block) {
			uint64_t i = index / N;
        	uint64_t k = index % N - i;
			calculate_dot_product(M, i, k, N);
			if (block) {
				bar.arrive_and_wait();
				if (PRINT_MESSAGES) std::printf("Computed index %lu and waited\n", index);
			}
			else {
				if (PRINT_MESSAGES) std::printf("Computed index %lu without waiting\n", index);
			}
		};

		auto wait = [&] () {
			bar.arrive_and_wait();
			if (PRINT_MESSAGES) std::printf("Waited\n");
		};

		ThreadPool TP(T);
		for (uint64_t k=0; k<N; ++k) { // for each upper diagonal
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
	int min                   = DEFAULT_MIN;
	int max                   = DEFAULT_MAX;
	uint64_t N                = DEFAULT_DIM;
	uint32_t T                = DEFAULT_NTHREADS;
	std::string log_file_name = DEFAULT_LOG_FILE;

	// Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 4 && argc != 5) {
        std::printf("use: %s N [min max] [T]\n", argv[0]);
        std::printf("     N size of the square matrix\n");
        std::printf("     min waiting time (us)\n");
        std::printf("     max waiting time (us)\n");
        std::printf("     T number of threads\n");
        return -1;
    }

	// Set the parameters (if any)
    if (argc > 1) {
        N = std::stol(argv[1]);
        if (argc > 2) {
            min = std::stol(argv[2]);
            max = std::stol(argv[3]);
            if (argc > 4) {
                T = std::stol(argv[4]);
            }
        }
    }

	// allocate the matrix
	std::vector<double> M(N*N, 0.0);

	uint64_t expected_totaltime=0;
	// init function
	auto init=[&]() {
		for(uint64_t k=0; k<N; ++k) {
			for(uint64_t i=0; i<(N-k); ++i) {
				if (k == 0) { // major diagonal
                    M[i * N + i] = (i + 1) / static_cast<double>(N);
                }
                else {
                    int t = random(min, max);
                    M[i * N + (i+k)] = t;
                    expected_totaltime += t;
                }
			}
		}
	};
	
	init();

	// sequential
	double sequential_time=-1;
	if (PRINT_MESSAGES) std::printf("------ Sequential Execution ------\n");
	TIMERSTART(wavefront_sequential);
	wavefront_sequential(M, N); 
	TIMERSTOP(wavefront_sequential, sequential_time);
	
	// parallel static
	double parallel_static_time=-1;
	if (PRINT_MESSAGES) std::printf("------ Parallel Static Execution ------\n");
	TIMERSTART(wavefront_parallel_static);
	wavefront_parallel_static(M, N, T); 
	TIMERSTOP(wavefront_parallel_static, parallel_static_time);

	// parallel dynamic
	double parallel_dynamic_time=-1;
	if (PRINT_MESSAGES) std::printf("------ Parallel Dynamic Execution ------\n");
	TIMERSTART(wavefront_parallel_dynamic);
	wavefront_parallel_dynamic(M, N, T); 
	TIMERSTOP(wavefront_parallel_dynamic, parallel_dynamic_time);


	// write the execution times to a file
	std::ofstream file;
	file.open(log_file_name, std::ios_base::app);
	file << N << "," << T << "," << min << "," << max << ","
		<< expected_totaltime/1000000 << "," << sequential_time << ","
		<< parallel_dynamic_time << "," << parallel_static_time << "\n";
	file.close();

	return 0;
}
