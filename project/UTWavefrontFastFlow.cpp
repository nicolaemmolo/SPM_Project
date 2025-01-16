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
#include <fstream>
#include <numeric>
#include <iomanip>
#include <hpc_helpers.hpp>
#include <random>


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
#define DEFAULT_MODE "ps"   // Default execution mode
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

/*
 * Class for the Wavefront computation using FastFlow
 */
class WavefrontFastFlow {
    // Emitter node: emits the diagonal index
    struct Emitter : ff_node_t<int> {
        uint64_t N;
        uint64_t current_diag;

        Emitter(uint64_t N) : N(N), current_diag(1) {}

        int* svc(int*) {
            if (current_diag >= N) return EOS; // Fine del lavoro
            auto diag = new int(current_diag); // Invio il numero della diagonale
            ++current_diag;
            return diag;
        }
    };

    // Worker node: computes the elements of the diagonal
    struct Worker : ff_node_t<int, std::pair<int, std::vector<double>>> {
        uint64_t N;
        std::vector<double>& M;

        Worker(uint64_t N, std::vector<double>& M) : N(N), M(M) {}

        std::pair<int, std::vector<double>>* svc(int* diag_ptr) {
            int diag = *diag_ptr;
            delete diag_ptr;

            uint64_t start = 0;
            uint64_t end = N - diag;
            std::vector<double> results(end - start);

            // Calcolo degli elementi della diagonale
            for (uint64_t i = start; i < end; ++i) {
                double result = 0.0;
                for (uint64_t j = 0; j < diag; ++j) {
                    result += M[INDEX(i, j, N)] * M[INDEX(diag + i - j, j, N)];
                }
                results[i - start] = std::cbrt(result); // Calcolo del risultato
            }
            return new std::pair<int, std::vector<double>>(diag, results);
        }
    };

    // Collector node: collects the results and updates the matrix
    struct Collector : ff_node_t<std::pair<int, std::vector<double>>> {
        uint64_t N;
        std::vector<double>& M;

        Collector(uint64_t N, std::vector<double>& M) : N(N), M(M) {}

        int* svc(std::pair<int, std::vector<double>>* result) {
            int diag = result->first;
            const std::vector<double>& values = result->second;

            for (size_t i = 0; i < values.size(); ++i) {
                M[INDEX(i, diag, N)] = values[i];
            }

            delete result;
            return GO_ON;
        }
    };
};


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
    Emitter emitter(N);
    Worker worker(N, M);
    Collector collector(N, M);

    ff_Farm<> farm;
    farm.add_emitter(emitter);
    farm.add_workers([&]() {
        std::vector<ff_node*> workers;
        for (uint32_t i = 0; i < T; ++i) {
            workers.push_back(new Worker(N, M));
        }
        return workers;
    }());
    farm.add_collector(collector);

    if (farm.run_and_wait_end() < 0) {
        std::cerr << "error during the execution (static)" << std::endl;
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
    uint32_t T                = DEFAULT_NTHREADS;
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
    file << N << "," << T << "," << mode << "," << execution_time << "\n";
    file.close();

    return 0;
}
