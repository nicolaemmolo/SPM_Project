#!/bin/bash
#SBATCH --job-name=wavefront_mpi    # Name of the job
#SBATCH --output=output.txt         # Name of output file
#SBATCH --error=error.txt           # Name of error file
#SBATCH --nodes=8                   # Number of computational nodes to be used
#SBATCH --ntasks-per-node=1         # Number of tasks per node
#SBATCH --time=01:59:00             # Maximum time requested for job execution (format HH:MM:SS)


# Parametri di test
#DIMENSIONS_LIST=(128 256 512 1024 2048 4096)
#NODES_LIST=(1 2 4 6 8)
#THREADS_LIST=(1 2 4 8 16 32)
#REPETITIONS=10

DIMENSIONS_LIST=(256 512)
NODES_LIST=(1 4)
THREADS_LIST=(1 2 4 8 16 32)
REPETITIONS=3

mpi_execution() { 
    for nodes in $($NODES_LIST); do
        for N in $($DIMENSIONS_LIST); do
            echo "N=$N results_UTWMPI.csv nodes=$nodes"
            for rep in $(seq 1 $REPETITIONS); do
                mpirun -n $nodes ./UTWMPI $N results_UTWMPI.csv $nodes
            done
        done
    done
}

mpi_omp_static_execution() {
    for nodes in $($NODES_LIST); do
        for N in $($DIMENSIONS_LIST); do
            for T in $($THREADS_LIST); do
                echo "N=$N T=$T ps results_UTWMPI_OMP_static.csv nodes=$nodes"
                for rep in $(seq 1 $REPETITIONS); do
                    mpirun -n $nodes ./UTWMPI_OMP $N $T ps results_UTWMPI_OMP_static.csv $nodes
                done
            done
        done
    done
}

mpi_omp_dynamic_execution() {
    for nodes in $($NODES_LIST); do
        for N in $($DIMENSIONS_LIST); do
            for T in $($THREADS_LIST); do
                echo "N=$N T=$T pd results_UTWMPI_OMP_dynamic.csv nodes=$nodes"
                for rep in $(seq 1 $REPETITIONS); do
                    mpirun -n $nodes ./UTWMPI_OMP $N $T pd results_UTWMPI_OMP_dynamic.csv $nodes
                done
            done
        done
    done
}


# MPI execution
mpi_execution

# MPI + OpenMP execution
#mpi_omp_static_execution
#mpi_omp_dynamic_execution