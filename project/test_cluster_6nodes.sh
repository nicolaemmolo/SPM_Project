#!/bin/bash
#SBATCH --job-name=wavefront_mpi    # Name of the job
#SBATCH --output=output6.txt        # Name of output file
#SBATCH --error=error6.txt          # Name of error file
#SBATCH --nodes=6                   # Number of computational nodes to be used
#SBATCH --time=01:59:00             # Maximum time requested for job execution (format HH:MM:SS)


# Test parameters
REPETITIONS=3


mpi_execution() {
    for N in 128 256 512 1024 2048 4096; do
        for rep in $(seq 1 $REPETITIONS); do
            mpirun ./UTWMPI $N results_UTWMPI_6.csv 6
        done
    done
}


# MPI execution
mpi_execution