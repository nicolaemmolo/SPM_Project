#!/bin/bash

# Parametri di test
NODES_LIST=(1 2 4 6 8)
REPETITIONS=10

# Funzione per sottomettere un lavoro SLURM
submit_job() {
    local NODES=$1
    local EXECUTABLE=$2
    local MODE=$3
    local OUTPUT_FILE=$4

    cat <<EOF | sbatch
#!/bin/bash
#SBATCH --job-name=${EXECUTABLE}_${MODE}_N${NODES}         # Nome del lavoro
#SBATCH --nodes=${NODES}                                 # Numero di nodi
#SBATCH --ntasks-per-node=8                              # Numero di task per nodo
#SBATCH --time=02:00:00                                 # Tempo massimo di esecuzione (HH:MM:SS)
#SBATCH --output=${OUTPUT_FILE}_%j.out                   # File di output
#SBATCH --error=${OUTPUT_FILE}_%j.err                    # File di errori

module load mpi                                        # Carica il modulo MPI
mpirun -n $(($NODES * 8)) ./${EXECUTABLE} 1000 ${MODE} ${OUTPUT_FILE}.csv   # Esegui il programma MPI
EOF
}

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
done