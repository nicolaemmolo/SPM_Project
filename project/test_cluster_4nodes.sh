#!/bin/bash
#SBATCH --job-name=wavefront_mpi    # Name of the job
#SBATCH --output=output4.txt        # Name of output file
#SBATCH --error=error4.txt          # Name of error file
#SBATCH --nodes=4                   # Number of computational nodes to be used
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --time=01:59:00             # Maximum time requested for job execution (format HH:MM:SS)


# Test parameters
REPETITIONS=3

mpi_execution() {
    for N in 128 256 512 1024 2048 4096; do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun ./UTWMPI $N results_UTWMPI_4.csv 4
        done
    done
}

mpi_omp_static_execution() {
    for N in 128 256 512 1024 2048 4096; do
        for T in 1 2 4 8 16 32; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun ./UTWMPI_OMP $N $T ps results_UTWMPI_OMP_static_4.csv 4
            done
        done
    done
}

mpi_omp_dynamic_execution() {
    for N in 128 256 512 1024 2048 4096; do
        for T in 1 2 4 8 16 32; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun ./UTWMPI_OMP $N $T pd results_UTWMPI_OMP_dynamic_4.csv 4
            done
        done
    done
}


# MPI execution
#mpi_execution

# MPI + OpenMP execution
mpi_omp_static_execution
#mpi_omp_dynamic_execution