//
// Sequential and Parallel code of the first SPM Assignment
//		with Wavefront computation
//
// compile:
// g++ -std=c++20 -O3 -march=native -Iinclude UTWavefrontMPI_OMP.cpp -o UTWMPIOMP
//

#include <iostream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <mpi.h>
#include <omp.h>

#define DEFAULT_DIM 3         // default size of the matrix (NxN)
#define DEFAULT_NTHREADS 2    // default number of OpenMP threads
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

// Main wavefront function with MPI and OpenMP
void wavefront_mpi_omp(std::vector<std::vector<double>> &M, const uint64_t &N, const int &rank, const int &size, const int &nthreads) {

    for (uint64_t k=1; k<N; ++k) { // for each upper diagonal
        #pragma omp parallel for num_threads(nthreads) schedule(static)
        for (uint64_t m=rank; m<(N-k); m+=size) { // distribute work among processes
            compute_diagonal_element(M, N, m, k);
        }

        // Synchronize all processes at the end of each diagonal step
        MPI_Barrier(MPI_COMM_WORLD);

        // Broadcast the matrix to ensure all processes have updated data
        for (uint64_t i = 0; i < N; ++i) {
            MPI_Bcast(&M[i][0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    uint64_t N                = DEFAULT_DIM;
    int nthreads              = DEFAULT_NTHREADS;
    std::string log_file_name = DEFAULT_LOG_FILE;

    // Verify the correct number of args
    if (argc != 1 && argc != 2 && argc != 3) {
        if (rank == 0) {
            std::printf("use: %s [N] [nthreads]\n", argv[0]);
            std::printf("     N        : size of the square matrix\n");
            std::printf("     nthreads : number of OpenMP threads per process\n");
        }
        MPI_Finalize();
        return -1;
    }

    // Set the parameters (if any)
    if (argc > 1) {
        N = std::stol(argv[1]);
        if (argc > 2) {
            nthreads = std::stoi(argv[2]);
        }
    }

    // allocate the matrix
    std::vector<std::vector<double>> M(N, std::vector<double>(N, 0.0));

    // init function
    auto init = [&]() {
        for (uint64_t m = 0; m < N; ++m) {
            M[m][m] = static_cast<double>(m + 1) / static_cast<double>(N);
        }
    };
    
    init();

    double execution_time = -1;

    if (rank == 0) std::printf("------ MPI + OpenMP Execution ------\n");

    double start_time = MPI_Wtime();
    wavefront_mpi_omp(M, N, rank, size, nthreads);
    double end_time = MPI_Wtime();

    execution_time = end_time - start_time;

    if (rank == 0 && PRINT_MATRIX) print_matrix(M, N);

    if (rank == 0) {
        // write the execution times to a file
        std::ofstream file;
        file.open(log_file_name, std::ios_base::app);
        file << N << "," << size << "," << nthreads << "," << "mpi_omp" << "," << execution_time << "\n";
        file.close();
    }

    MPI_Finalize();
    return 0;
}
