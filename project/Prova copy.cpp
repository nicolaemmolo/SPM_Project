#include <iostream>         // For input/output stream
#include <vector>           // For using the std::vector container
#include <cmath>            // For mathematical functions, here specifically for std::cbrt
#include <ff/ff.hpp>        // For FastFlow framework
#include <ff/farm.hpp>      // For FastFlow farm
#include <hpc_helpers.hpp>  // This file include TIMER macros
#include <barrier>          // C++20 barrier for thread synchronization

using namespace ff;

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


        // Process each upper diagonal
        for (uint64_t k = 1; k < N; ++k) {
            // Recompute variable for k-th diagonal
            uint64_t chunk_size = (N-k) / T;
            uint64_t remainder = (N-k) % T;

            // Recompute the interval for each worker
            start = task->rank * chunk_size + (task->rank < remainder ? task->rank : remainder);
            end = (task->tank + 1) * chunk_size + (task->rank < remainder ? (task->rank + 1) : remainder);

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
void wavefront(std::vector<double> &M, const uint64_t &N, const uint64_t &T) {
    std::vector<Task> tasks;

    std::barrier barrier(T);

    // Compute interval values for the first upper diagonal
    uint64_t chunk_size = (N-1) / T; // Compute chunk size
    uint64_t remainder = (N-1) % T; // If N-1 is not divisible by T then there will be a remainder
    uint64_t start = 0;

    // Create as many tasks as workers to use
    for (uint64_t int t = 0; t < T; ++w) {
        uint64_t end = start + chunk_size + (t < remainder ? 1 : 0);
        tasks.push_back(Task{start, end, N, T, t, &M});
        start = end;
    }

    // Pass tasks to emitter
    Emitter emitter(tasks);
    std::vector<std::unique_ptr<ff_node>> workers;
    for (uint64_t int i = 0; i < T; ++i) {
        workers.push_back(make_unique<Worker>(barrier));
    }

    // Create MAP
    ff_Farm<Task> map(std::move(workers), emitter);
    map.remove_collector();    // Remove collector as we don't need to collect results
    map.set_scheduling_ondemand(); // Set scheduling policy

    // run Map
    if (map.run_and_wait_end() < 0) {
        error("running farm");
    }
}



int main(int argc, char *argv[]) {
    uint64_t N = 512; // Default size of the matrix (NxN)
    unsigned int num_workers = ff_numCores(); // Default number of workers
    
    if (argc < 2 || argc > 3) {
        std::printf("use: %s N [num_workers]\n", argv[0]);
        std::printf("     N size of the square matrix\n");
        std::printf("     num_workers number of workers to use\n");
        return -1;
    }
    if (argc > 1) {
        N = std::stol(argv[1]);
    }
    if (argc > 2) {
        num_workers = std::stoul(argv[2]);
    }

    // Allocate and initialize the matrix
    std::vector<std::vector<double>> M(N, std::vector<double>(N, 0));
    auto init = [&]() {
        for (uint64_t i = 0; i < N; ++i) {
            M[i][i] = static_cast<double>(i + 1) / N;
        }
    };
    init();

    std::printf("Matrix initialized.\n");

    // Measure the time taken by the wavefront computation
    TIMERSTART(wavefront);
    wavefront(M, N, num_workers);
    TIMERSTOP(wavefront);

    // Print the value in the upper right corner of the matrix
    std::cout << "Value in the upper right corner of the matrix: " << M[0][N-1] << std::endl;
    return 0;
}