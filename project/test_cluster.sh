#!/bin/bash
#SBATCH --job-name=wavefront_mpi   # Nome del job
#SBATCH --output=output.txt        # Nome del file di output
#SBATCH --error=error.txt          # Nome del file di errore
#SBATCH --ntasks=8                 # Numero totale di task MPI
#SBATCH --nodes=2                  # Numero di nodi richiesti
#SBATCH --time=00:10:00            # Tempo massimo di esecuzione (hh:mm:ss)
#SBATCH --partition=normal         # Partizione o coda da usare (sostituisci con la partizione corretta del tuo cluster)

#module load mpi                    # Carica il modulo MPI, se necessario

mpirun ./UTWMPI 3 wavefront_results_MPI.csv 2


# Parametri di test
NODES_LIST=(1 2 4 6 8)
REPETITIONS=10
A = '
# Test per UTWavefrontMPI
for NODES in "${NODES_LIST[@]}"; do
    for ((rep=1; rep<=REPETITIONS; rep++)); do
        echo "Testing UTWavefrontMPI with ${NODES} nodes, repetition ${rep}"
        submit_job ${NODES} UTWavefrontMPI s results_UTWavefrontMPI
    done
done

# Test per UTWavefrontMPI_FF
for NODES in "${NODES_LIST[@]}"; do
    for ((rep=1; rep<=REPETITIONS; rep++)); do
        echo "Testing UTWavefrontMPI_FF with ${NODES} nodes, repetition ${rep}"
        submit_job ${NODES} UTWavefrontMPI_FF ps results_UTWavefrontMPI_FF_static
        submit_job ${NODES} UTWavefrontMPI_FF pd results_UTWavefrontMPI_FF_dynamic
    done
done

# Test per UTWavefrontMPI_OMP
for NODES in "${NODES_LIST[@]}"; do
    for ((rep=1; rep<=REPETITIONS; rep++)); do
        echo "Testing UTWavefrontMPI_OMP with ${NODES} nodes, repetition ${rep}"
        submit_job ${NODES} UTWavefrontMPI_OMP ps results_UTWavefrontMPI_OMP_static
        submit_job ${NODES} UTWavefrontMPI_OMP pd results_UTWavefrontMPI_OMP_dynamic
    done
done'