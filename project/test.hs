#!/bin/bash

THREADS_STEP=2
MAX_THREADS=40
REPETITIONS=10


# Testing Sequential Execution
seq_execution() {
    for N in 128 256 512 1024 2048 4096; do
        echo "N=$N T=1 s results_UTW_seq.csv"
        for rep in $(seq 1 $REPETITIONS); do
            ./UTW $N 1 s results_UTW.csv
        done
    done
}

# Testing Parallel Static Execution
par_static_execution() {
    for T in $(seq 4 $THREADS_STEP $MAX_THREADS); do
        for N in 128 256 512 1024 2048; do
            echo "N=$N T=$T ps results_UTW_static.csv"
            for rep in $(seq 1 $REPETITIONS); do
                ./UTW $N $T ps results_UTW_static.csv
            done
        done
    done
}

# Testing Parallel Dynamic Execution
par_dynamic_execution() {
    for T in $(seq 4 $THREADS_STEP $MAX_THREADS); do
        for N in 128 256 512 1024 2048; do
            echo "N=$N T=$T pd results_UTW_dynamic.csv"
            for rep in $(seq 1 $REPETITIONS); do
                ./UTW $N $T pd results_UTW_dynamic.csv
            done
        done
    done
}

# Execution
#seq_execution
par_static_execution
#par_dynamic_execution
