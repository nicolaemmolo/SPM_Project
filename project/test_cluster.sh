#!/bin/bash
#SBATCH --job-name=wavefront_mpi    # Name of the job
#SBATCH --output=output.txt         # Name of output file
#SBATCH --error=error.txt           # Name of error file
#SBATCH --nodes=8                   # Number of computational nodes to be used
#SBATCH --time=01:59:00             # Maximum time requested for job execution (format HH:MM:SS)


# Test parameters
REPETITIONS=10


# Testing MPI Execution (with 1 node)
mpi_execution_1() {
    for N in 128 256 512 1024 2048 4096; do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun -n 1 ./UTWMPI $N results_UTWMPI_1.csv 1
        done
    done
}

# Testing MPI Execution (with 2 nodes)
mpi_execution_2() {
    for N in 128 256 512 1024 2048 4096; do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun -n 2 ./UTWMPI $N results_UTWMPI_2.csv 2
        done
    done
}

# Testing MPI Execution (with 4 nodes)
mpi_execution_4() {
    for N in 128 256 512 1024 2048 4096; do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun -n 4 ./UTWMPI $N results_UTWMPI_4.csv 4
        done
    done
}

# Testing MPI Execution (with 6 nodes)
mpi_execution_6() {
    for N in 128 256 512 1024 2048 4096; do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun -n 6 ./UTWMPI $N results_UTWMPI_6.csv 6
        done
    done
}

# Testing MPI Execution (with 8 nodes)
mpi_execution_8() {
    for N in 128 256 512 1024 2048 4096; do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun -n 8 ./UTWMPI $N results_UTWMPI_8.csv 8
        done
    done
}

# Testing MPI+OMP Execution (with 1 node)
mpi_omp_execution_1() {
    for T in 1 2 4 8 16 32; do
        for N in 128 256 512 1024 2048 4096; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun -n 1 ./UTWMPIOMP $N results_UTWMPIOMP_1.csv 1 $T
            done
        done
    done
}

# Testing MPI+OMP Execution (with 2 nodes)
mpi_omp_execution_2() {
    for T in 1 2 4 8 16 32; do
        for N in 128 256 512 1024 2048 4096; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun -n 2 ./UTWMPIOMP $N results_UTWMPIOMP_2.csv 2 $T
            done
        done
    done
}

# Testing MPI+OMP Execution (with 4 nodes)
mpi_omp_execution_4() {
    for T in 1 2 4 8 16 32; do
        for N in 128 256 512 1024 2048 4096; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun -n 4 ./UTWMPIOMP $N results_UTWMPIOMP_4.csv 4 $T
            done
        done
    done
}

# Testing MPI+OMP Execution (with 6 nodes)
mpi_omp_execution_6() {
    for T in 1 2 4 8 16 32; do
        for N in 128 256 512 1024 2048 4096; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun -n 6 ./UTWMPIOMP $N results_UTWMPIOMP_6.csv 6 $T
            done
        done
    done
}

# Testing MPI+OMP Execution (with 8 nodes)
mpi_omp_execution_8() {
    for T in 1 2 4 8 16 32; do
        for N in 128 256 512 1024 2048 4096; do
            for rep in $(seq 1 $REPETITIONS); do
                mpirun -n 8 ./UTWMPIOMP $N results_UTWMPIOMP_8.csv 8 $T
            done
        done
    done
}


# MPI execution
#mpi_execution_1
#mpi_execution_2
#mpi_execution_4
#mpi_execution_6
#mpi_execution_8

# MPI+OMP execution
#mpi_omp_execution_1
#mpi_omp_execution_2
#mpi_omp_execution_4
#mpi_omp_execution_6
#mpi_omp_execution_8