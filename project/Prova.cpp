#include <vector>
#include <random>
#include <ff/ff.hpp>

#include <cmath>
#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <filesystem>

using namespace ff;

// mutex only used if debug mode is used
#include <mutex>
std::mutex mtx;



// MACROS

#define DEFAULT_INITIALIZATION 0.0 // default matrix value inzialization (double)s
#define DEFAULT_M_SIZE 10 // default size of the square matrix (NxN)
#define DEFAULT_GRANULARITY 3 // default size of the task (NxN)
#define DEFAULT_W 2 // default number of workers
#define DEFAULT_LOG_NAME "./log/fastflow_stream_sm.csv" // default log file name
#define DEFAULT_MACHINE 0 // default machine identifier



// Define types for task source and task sink, the elements that source-sink and worker will send to each other
// task_source contains information to process a part of the matrix
using task_source = std::tuple<TriangularMatrix<double>*, int, int, int, int>;
// task_sink is used to indicate completion of a processing task
using task_sink = std::tuple<int, int>;

// Class to handle the source-sink operations, sending task to the workers and updating the "global" matrix
struct SourceSinkSingleMatrix: ff_monode_t<task_sink, task_source> {

    // Initial value used to populate the matrix
    const double init_value = DEFAULT_INITIALIZATION;
    
    // Matrix where to store results and matrix where to keep track of processed elements
    TriangularMatrix<double> row_tm;
    TriangularMatrix<bool> presence_um;
    
    // Dimensions and granularity of the matrix
    const int true_matrix_dimension;
    const int pres_matrix_dimension;
    // Used to determine the task dimension and to divide the global task
    const int granularity;

    // Constructor to initialize the SourceSinkSingleMatrix
    SourceSinkSingleMatrix(const int matrix_dimension, const int granularity, double* diagonal):
        // Initialisation real matrix, element matrix
        row_tm(matrix_dimension, init_value),
        // Initialisation virtual matrix, task matrix
        presence_um(std::ceil(static_cast<double>(matrix_dimension - 1)/granularity), 0),
        // Dimension (number of row or column) of the real element matrix
        true_matrix_dimension(matrix_dimension),
        // Dimension (number of row or column) of the presence matrix, virtual matrix, task matrix
        pres_matrix_dimension(std::ceil(static_cast<double>(matrix_dimension - 1)/granularity)),
        granularity(granularity)
        {
            // Output debug information abount sink initialization, level 3
            DEBUG_STDOUT( 3, "Sink node initialized with dim: " << matrix_dimension )
            
            // Initialize the diagonal elements of the matrix
            for (int k = 0; k < true_matrix_dimension; k++){
                row_tm.element(k,k) = diagonal[k];
            }
            // Free the memory allocated for the diagonal array
            delete[] diagonal;
        }

    // Called when the node is initializated, send initial tasks to generate the "sub-task" diagonal
    int svc_init(){
        DEBUG_STDOUT( 3, "sink svc_init()" )

        // Send tasks for each row of the matrix
        for (int k = 0; k < pres_matrix_dimension; k++){
            ff_send_out(new task_source(&row_tm, k, k, granularity, 1));
        }
        return 0; 
    }

    // Process received (completed) tasks
    task_source* svc(task_sink* task) {

        // Handle null tasks sended by fastflow
        if (task == nullptr) {
            DEBUG_STDOUT(5, "SOURCE: null "<< task)
            return this->GO_ON;
        }

        // tuple unpacking
        int i = std::get<0>(*task);
        int j = std::get<1>(*task);
        DEBUG_STDOUT(5, "SOURCE: " << i << " - " << j)

        // Mark the received task as processed
        presence_um.element(i,j) = 1;

        // tasks are sent if they can be calculated, depending on the presence of the two adjacent tasks.
        if (i-1 >= 0  && j-1 >= 0 && presence_um.element(i-1, j-1))
            ff_send_out(new task_source(&row_tm, i-1 , j, granularity, 0));
        if (i+1 < pres_matrix_dimension  && j+1 < pres_matrix_dimension && presence_um.element(i+1, j+1))
            ff_send_out(new task_source(&row_tm, i , j+1, granularity, 0));

        // Broadcast end-of-stream signal if is reached the end of the matrix 
        if (i == 0 && j == pres_matrix_dimension - 1){
            broadcast_task(EOS);
            return EOS;
        }
        else return this->GO_ON;
    }

    // Printing debug information at the termination
    void svc_end(){
        DEBUG_STDOUT(3, "sink svc_end()")
        DEBUG_STDOUT(5, "Final Matrix row:" << std::endl << row_tm)
        DEBUG_STDOUT(2, "Final value (matrix dim " << true_matrix_dimension << "): " << row_tm.element(0, true_matrix_dimension-1))
    } 

    // Printing debug information
    void eosnotify(ssize_t ch){
        DEBUG_STDOUT(3, "sink eosnotify()")
    } 
};

// Class that implements the worker nodes
struct WorkerSingleMatrix: ff_node_t<task_source, task_sink> {
    
    // Compute one element of the wavefront computation that comprise the vectorial sum 
    // and the cubic root
    double compute(TriangularMatrix<double>* tm, int i, int j){
        double partial_sum = 0;
        for (int k=0; k < (j-i); k++){
            // Compute partial sum for the row and column adjacent vector
            partial_sum += tm->element(i, i + k) * tm->element(j - k, j);
        }
        // Return the cubic root of the partial sum
        return std::cbrt(partial_sum);
    }
    
    // Process incoming tasks sent by the source-sink
    task_sink* svc(task_source* task) {

        // tuple unpacking
        TriangularMatrix<double>* tm = std::get<0>(*task);
        int i = std::get<1>(*task);
        int j = std::get<2>(*task);
        int granularity = std::get<3>(*task);
        int first_it = std::get<4>(*task);

        // Print debug information if DEBUG > 5
        DEBUG_STDOUT(5, "WORKER: " << i << " - " << j << " - " << first_it)

        // Variables for the coordinates relative to the element in the matrix, not the task matrix
        int abs_i;
        int abs_j;

        // If this is NOT the first iteration, calculate the lower triangular part of the task,
        // if this is the firts iteration, the lower triangular matrix is not needed and not computable
        if (!first_it) {
            // Iterate over diagonal
            for (int diag = 0; diag < granularity - 1; diag++){      
                // Iterate over diagonal element
                for (int k=diag; k>=0; k--){
                    // Compute the absolute (element matrix) coordinates
                    abs_i = i*granularity + (granularity - k - 1);
                    abs_j = j*granularity + 1 + (diag - k);

                    // Break if out of bounds
                    if (abs_i >= tm->matrix_dimension || abs_j >= tm->matrix_dimension) break;

                    // Update element matrix with computed value
                    tm->element(abs_i, abs_j) = compute(tm, abs_i, abs_j);
                }
            }
        }

        // Compute upper triangular section of the task for all iteration 
        // (this use the lower part if is not the first iteration)
        // Iterate over diagonal
        for (int diag = 0; diag < granularity; diag++){
            // Iterate over diagonal element
            for (int k=0; k < (granularity - diag); k++){
                // Compute the absolute (element matrix) coordinates
                abs_i = i*granularity + (k);
                abs_j = j*granularity + 1 + (diag + k);

                // Break if out of bounds
                if (abs_i >= tm->matrix_dimension || abs_j >= tm->matrix_dimension) 
                    break;

                // Update element matrix with computed value
                tm->element(abs_i, abs_j) = compute(tm, abs_i, abs_j);
            }
        }

        // Send completed task to source-sink
        ff_send_out(new task_sink(i, j));
        delete task; // Clean up task memory
        return GO_ON; // Continue processing; do not terminate
    }

    // Printing debug information
    int svc_init(){
        DEBUG_STDOUT(3, "worker svc_init()")
        return 0; 
    }
    
    // Printing debug information
    void svc_end(){
        DEBUG_STDOUT(3, "worker svc_end()")
    } 
    
    // Printing debug information
    void eosnotify(ssize_t ch){
        DEBUG_STDOUT(3, "worker eosnotify()")
    } 
};





// MAIN

int main(int argc, char *argv[]) {   

    // Set default values for matrix size, granularity, number of workers, machine identifier, and log file name
    int matrix_size = DEFAULT_M_SIZE;
    int granularity = DEFAULT_GRANULARITY;
    int num_workers = DEFAULT_W;
    int machine_identifier = DEFAULT_MACHINE;
    std::string log_file_name = DEFAULT_LOG_NAME;


    // Record initialization time
    double init_time;
    TIMERSTART(init)
    // Initialize matrix diagonal values
    double* init = new double[matrix_size];
    for (int i = 0; i<matrix_size; i++){
        init[i] = (i+1.0)/matrix_size;
    }

    // Create the Source-Sink object with initial values
    SourceSinkSingleMatrix producer(matrix_size, granularity, init);

    // Configure the FastFlow farm with custon worker nodes
    ff_Farm<task_source, task_sink> farm(
                            [&]() {
                                std::vector<std::unique_ptr<ff_node> > W;
                                // Initialize worker nodes
                                for(auto i=0;i<num_workers;++i)
                                    W.push_back(make_unique<WorkerSingleMatrix>());
                                return W;
                            } ()
                            ,
                            producer);

    // Remove the default collector from the farm as it's not needed
    farm.remove_collector();
    // Wrap the farm to enable feedback from Workers to the Emitter
    farm.wrap_around();        

    // Stop the initialization timer and print elapsed time
    TIMERSTOP(init, init_time)
  
  

    // Record running time
    double running_time;
    TIMERSTART(run)
    // Run the FastFlow farm and wait for its completion
    if (farm.run_and_wait_end()<0) {
        std::cerr << "Error running farm" << std::endl;
        return -1;
    }
    // Stop the computation timer and print elapsed time
    TIMERSTOP(run, running_time)
    DEBUG_STDOUT(2,
                "Computation complete - elapsed time: " << running_time*1000
            )

    // Log the execution times and configuration parameters to the log file
    std::ofstream file;
    file.open(log_file_name, std::ios_base::app);
    file << matrix_size << "," << num_workers << "," << machine_identifier << "," << granularity << "," << DEBUG << ","
        << init_time*1000  << "," << running_time*1000 << std::endl;
    file.close();

    return 0;
}