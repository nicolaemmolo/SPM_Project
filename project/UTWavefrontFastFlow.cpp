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
#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>
#include <ff/farm.hpp>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <hpc_helpers.hpp>
#include <random>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <barrier>


using namespace ff;

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
#define DEFAULT_NTHREADS 2  // Default number of threads
#define DEFAULT_MODE "f"   // Default execution mode
#define DEFAULT_LOG_FILE "wavefront_results.csv"    // Default log file name

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


// Define a struct to represent a task that each worker will execute
struct Task {
    uint64_t start, end, N, T, rank;
    std::vector<double>* M;
};

// Define a worker class that inherits from ff_node_t and processes tasks
struct Worker : ff_node_t<Task> {
    Worker(std::barrier<> &barrier) : barrier(barrier) {}

    Task* svc(Task* task) {
        auto& M = *task->M; // Matrix on which computation takes
        auto N = task->N; // Size of the matrix
        auto T = task->T; // Number of threads
        auto start = task->start; // Start index of the diagonal
        auto end = task->end; // End index of the diagonal

        // Process each upper diagonal
        for (uint64_t k = 1; k < N; ++k) {
            // Recompute variable for k-th diagonal
            uint64_t chunk_size = (N-k) / T;
            uint64_t remainder = (N-k) % T;

            // Recompute the interval for each worker
            start = task->rank * chunk_size + (task->rank < remainder ? task->rank : remainder);
            end = (task->rank + 1) * chunk_size + (task->rank < remainder ? (task->rank + 1) : remainder);

            // Process elements in the k-th diagonal assigned to this worker
            for (uint64_t i = start; i < end; ++i) {
                compute_diagonal_element(M, N, i, k);
            }
            barrier.arrive_and_wait(); // Synchronize after processing each diagonal
        }
        delete task; // Clean up task after processing
        return GO_ON;
    }

    std::barrier<> &barrier;
};

// Define an emitter class that inherits from ff_monode_t and generates tasks
struct Emitter : ff_monode_t<Task> {
    Emitter(const std::vector<Task>& tasks) : tasks(tasks), task_index(0) {}

    Task* svc(Task*) {
        if (task_index >= tasks.size())
            return EOS;
        return new Task(tasks[task_index++]);
    }

    std::vector<Task> tasks;
    size_t task_index;
};

// Function to perform wavefront computation on matrix M of size N with num_workers
void wavefront_farm(std::vector<double> &M, const uint64_t &N, const uint64_t &T) {
    std::vector<Task> tasks;

    std::barrier barrier(T);

    // Compute interval values for the first upper diagonal
    uint64_t chunk_size = (N-1) / T; // Compute chunk size
    uint64_t remainder = (N-1) % T; // If N-1 is not divisible by T then there will be a remainder
    uint64_t start = 0;

    // Create as many tasks as workers to use
    for (uint64_t t = 0; t < T; ++t) {
        uint64_t end = start + chunk_size + (t < remainder ? 1 : 0);
        tasks.push_back(Task{start, end, N, T, t, &M});
        start = end;
    }

    // Pass tasks to emitter
    Emitter emitter(tasks);
    std::vector<std::unique_ptr<ff_node>> workers;
    for (uint64_t i = 0; i < T; ++i) {
        workers.push_back(make_unique<Worker>(barrier));
    }

    // Create Farm
    ff_Farm<Task> farm(std::move(workers), emitter);
    farm.remove_collector();    // Remove collector as we don't need to collect results
    farm.set_scheduling_ondemand(); // Set scheduling policy

    // Run Farm
    if (farm.run_and_wait_end() < 0) {
        error("running farm");
    }
}


/* Wavefront (parallel version with static scheduling using FastFlow)
 * @param M: matrix
 * @param N: size of the matrix
 * @param T: number of threads
 */
void wavefront_parallel_static_ff(std::vector<double> &M, const uint64_t &N, const uint32_t &T) {
    ParallelFor pf(T);
    
	for (uint64_t k=1; k<N; ++k) { // For each upper diagonal
        pf.parallel_for(0, N-k, 1, 0, [&](const long i) {
            compute_diagonal_element(M, N, i, k);
        });
    }
}


/* Wavefront (parallel version with dynamic scheduling using FastFlow)
 * @param M: matrix
 * @param N: size of the matrix
 * @param T: number of threads
 */
void wavefront_parallel_dynamic_ff(std::vector<double> &M, const uint64_t &N, const uint32_t &T) {
    ParallelFor pf(T);
    
	for (uint64_t k=1; k<N; ++k) { // For each upper diagonal
        pf.parallel_for(0, N-k, 1, 1, [&](const long i) {
            compute_diagonal_element(M, N, i, k);
        });
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


    // Allocate the matrix as a 1D array
    uint64_t total_elements = (N*(N+1))/2; // Total number of elements in the upper triangular matrix
    std::vector<double> M(total_elements, 0.0);

    // Init function (initialize the diagonal elements)
	auto init=[&]() {
		for (uint64_t i = 0; i < N; ++i) {
				M[i] = static_cast<double>(i+1) / static_cast<double>(N);
			}
	};
    
    init();

    double execution_time=-1;

    // Parallel farm
    if (mode == "f") {
        if (PRINT_MESSAGE) std::printf("------ Farm Execution ------\n");
        TIMERSTART(wavefront_farm);
        wavefront_farm(M, N, T);
        TIMERSTOP(wavefront_farm, execution_time);
    }

    // Parallel static
	if (mode == "ps") {
		if (PRINT_MESSAGE) std::printf("------ Parallel Static Execution ------\n");
		TIMERSTART(wavefront_parallel_static_ff);
		wavefront_parallel_static_ff(M, N, T); 
		TIMERSTOP(wavefront_parallel_static_ff, execution_time);
	}

	// Parallel dynamic
	if (mode == "pd") {
		if (PRINT_MESSAGE) std::printf("------ Parallel Dynamic Execution ------\n");
		TIMERSTART(wavefront_parallel_dynamic_ff);
		wavefront_parallel_dynamic_ff(M, N, T); 
		TIMERSTOP(wavefront_parallel_dynamic_ff, execution_time);
	}

	if (PRINT_MATRIX) print_matrix(M,N);
    if (PRINT_MATRIX) print_M(M,total_elements);
	if (PRINT_LAST_ELEMENT) print_last_element(M,total_elements);

    // Write the execution times to a file
    std::ofstream file;
    file.open(log_file_name, std::ios_base::app);
    file << N << "," << T << "," << mode << "," << execution_time << print_last_element(M,total_elements) << "\n";
    file.close();

    return 0;
}