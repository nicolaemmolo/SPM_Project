#!/bin/bash
#SBATCH --job-name=wavefront_mpi    # Name of the job
#SBATCH --output=output1.txt         # Name of output file
#SBATCH --error=error1.txt           # Name of error file
#SBATCH --nodes=1                   # Number of computational nodes to be used
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --time=00:10:00             # Maximum time requested for job execution (format HH:MM:SS)


# Parametri di test
DIMENSIONS_LIST=(128 256 512 1024 2048 4096)
THREADS_LIST=(1 2 4 8 16 32)
REPETITIONS=3

mpi_execution() {
    for N in $($DIMENSIONS_LIST); do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun ./UTWMPI $N results_UTWMPI_1.csv $nodes
        done
    done
}

mpi_omp_static_execution() {
    for N in $($DIMENSIONS_LIST); do
        for T in $($THREADS_LIST); do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun ./UTWMPI_OMP $N $T ps results_UTWMPI_OMP_static_1.csv $nodes
            done
        done
    done
}

mpi_omp_dynamic_execution() {
    for N in $($DIMENSIONS_LIST); do
        for T in $($THREADS_LIST); do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun ./UTWMPI_OMP $N $T pd results_UTWMPI_OMP_dynamic_1.csv $nodes
            done
        done
    done
}


# MPI execution
mpi_execution

# MPI + OpenMP execution
#mpi_omp_static_execution
#mpi_omp_dynamic_execution