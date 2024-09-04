#!/bin/bash
#SBATCH --job-name=wavefront_mpi   # Nome del job
#SBATCH --output=output.txt        # Nome del file di output
#SBATCH --error=error.txt          # Nome del file di errore
#SBATCH --ntasks=8                 # Numero totale di task MPI
#SBATCH --nodes=2                  # Numero di nodi richiesti
#SBATCH --time=00:10:00            # Tempo massimo di esecuzione (hh:mm:ss)
#SBATCH --partition=normal         # Partizione o coda da usare (sostituisci con la partizione corretta del tuo cluster)

mpirun -n 1 ./UTWMPI 128 wavefront_results_MPI.csv 2
mpirun -n 1 ./UTWMPIOMP 128 2 ps wavefront_results_MPI.csv 2

# Parametri di test
DIMENSIONS_LIST=(128 256 512 1024 2048 4096)
NODES_LIST=(1 2 4 6 8)
THREADS_LIST=(1 2 4 8 16 32)
REPETITIONS=10

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
#mpi_execution

# MPI + OpenMP execution
#mpi_omp_static_execution
#mpi_omp_dynamic_execution