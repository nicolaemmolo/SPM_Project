# SPM Project: Wavefront Computation Algorithm

Parallel and Distributed Systems project. Computer Science Master Degree, University of Pisa. A.Y 2023/2024


## Overview & Project Goal

This repository contains the final project for the Parallel and Distributed Systems course. The main objective was to design, implement, and evaluate various parallelization strategies for an algorithm following a **Wavefront Computation pattern** on a square matrix (NxN). Parallelization approaches using FastFlow, OpenMP and MPI were used, with tests performed on both a single multi-core machine and a cluster of multi-core machines. The aim of the project was to measure and compare the performance of each approach in terms of execution time, acceleration, efficiency and scalability.

**The Problem:**
In a wavefront pattern, the elements of a matrix are computed diagonally, starting from the top-left main diagonal and moving towards the top-right corner. 
* **Parallelism:** Elements lying on the *same* diagonal are strictly independent and can be computed concurrently.
* **Dependency:** A diagonal `k+1` cannot be started until the previous diagonal `k` is fully computed.
* **The Computation:** To move beyond a simple "time-wasting" assignment, the value of each element in this project is calculated as the cubic root of the dot product between its corresponding row and column vectors.

**The Goal:**
The project aims to overcome the strict synchronization barriers of this pattern by exploiting different parallel programming frameworks. The solutions are evaluated based on their **Execution Time, Speedup, Efficiency, and Scalability (Strong and Weak)** across both single multi-core machines and distributed clusters.

## Parallelization Strategies

The project explores multiple frameworks to handle the workload efficiently:

**1. Single Multi-Core Machine:**
* **C++ Native Threads:** Custom parallelization using `std::thread` and the C++20 `std::barrier` for strict diagonal synchronization (Static and Dynamic scheduling).
* **FastFlow:** Advanced task-based parallelism utilizing the `Farm` pattern and `parallel_for` features.
* **OpenMP:** Compiler-directive-based parallelization for quick and standardized multi-threading.

**2. Cluster of Multi-Core Machines:**
* **MPI (Message Passing Interface):** Distributed memory parallelization, splitting the matrix across multiple network nodes and syncing data via `MPI_Allgatherv`.
* **Hybrid MPI + OpenMP:** A mixed approach combining distributed memory (MPI for nodes) and shared memory (OpenMP for internal node threads).


## Directory Structure (inside `project/` folder)

### Source Code & Headers
* **`include/`**: Contains custom header files required for the project.
* **`UTWavefront.cpp`**: Baseline sequential implementation and native C++ Threads version.
* **`UTWavefrontFastFlow.cpp`**: Parallel implementation using the FastFlow library.
* **`UTWavefrontOpenMP.cpp`**: Parallel implementation using OpenMP directives.
* **`UTWavefrontMPI.cpp`**: Distributed parallel implementation using MPI.
* **`UTWavefrontMPIOMP.cpp`**: Hybrid distributed implementation using MPI + OpenMP.

### Automation & Scripts
* **`Makefile`**: Automates the compilation process for all versions (FastFlow, OpenMP, MPI).
* **`test_numa.sh`**: Bash script to automate performance testing on the local NUMA multi-core machine.
* **`test_cluster.sh`**: Bash script to automate testing across a different number of nodes on the cluster.

### Analysis & Results
* **`results/`**: Directory storing the raw output logs of the executed tests.
* **`plots_numa.ipynb`**: Jupyter Notebook for analyzing and visualizing (Speedup, Scalability) single-machine results.
* **`plots_cluster.ipynb`**: Jupyter Notebook for visualizing distributed cluster results.

### Executables (Generated via Make)
* **`UTW`**: Executable for Sequential and C++ Threads.
* **`UTWFF`**: Executable for FastFlow.
* **`UTWOMP`**: Executable for OpenMP.
* **`UTWMPI`**: Executable for MPI.
* **`UTWMPIOMP`**: Executable for MPI + OpenMP.
