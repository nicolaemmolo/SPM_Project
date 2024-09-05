#!/bin/bash
#SBATCH --job-name=wavefront_mpi    # Name of the job
#SBATCH --output=output11.txt        # Name of output file
#SBATCH --error=error11.txt          # Name of error file
#SBATCH --nodes=1                   # Number of computational nodes to be used
#SBATCH --time=01:59:00             # Maximum time requested for job execution (format HH:MM:SS)


# Test parameters
REPETITIONS=3

#mpirun ./UTWMPI_OMP $N $T ps results_UTWMPI_OMP_static_1.csv 1

mpi_execution() {
    for N in 128 256 512 1024 2048 4096; do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun ./UTWMPI $N results_UTWMPI_1.csv 1
        done
    done
}

mpi_omp_static_execution() {
    for N in 128 256 512 1024 2048 4096; do
        for T in 1 2 4 8 16 32; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun ./UTWMPIOMP $N $T ps results_UTWMPI_OMP_static_1.csv 1
            done
        done
    done
}

mpi_omp_dynamic_execution() {
    for N in 128 256 512 1024 2048 4096; do
        for T in 1 2 4 8 16 32; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun ./UTWMPIOMP $N $T pd results_UTWMPI_OMP_dynamic_1.csv 1
            done
        done
    done
}


# MPI execution
#mpi_execution

# MPI + OpenMP execution
#mpi_omp_static_execution
#mpi_omp_dynamic_execution