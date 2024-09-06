# SPM Project: Wavefront Computation Algorithm

## Overview

This project implements various parallelisation techniques for a wavefront algorithm applied to square matrices. Parallelisation approaches using FastFlow, OpenMP and MPI were used, with tests performed on both a single multi-core machine and a cluster of multi-core machines. The aim of the project was to measure and compare the performance of each approach in terms of execution time, acceleration, efficiency and scalability.


## Directory Structure

* **include/**: Contains header files for the project.
* **results/**: Directory where the results of executed tests are stored.
* **Makefile**: File for the automatic compilation of the sources. Supports the creation of executables for the different implementations (FastFlow, OpenMP, MPI).
* **plots_cluster.ipynb**: Jupyter notebook for displaying the results of tests performed on the cluster.
* **plots_numa.ipynb**: Jupyter notebook for displaying results of tests run on a single NUMA machine.
* **test_cluster_*.sh**: Script for running automated tests on the cluster with a different number of nodes (1, 2, 4, 6, 8 nodes).
* **test_numa.sh**: Script to run tests on the local NUMA multi-core machine.
* **UTWavefront.cpp**: Basic sequential implementation of the UTWavefront algorithm.
* **UTWavefrontFastFlow.cpp**: Implementation of the UTWavefront algorithm using the FastFlow library for parallelisation.
* **UTWavefrontOpenMP.cpp**: Implementation of the UTWavefront algorithm using OpenMP for parallelisation.
* **UTWavefrontMPI.cpp**: Implementation of the UTWavefront algorithm using MPI for parallelisation on a cluster of machines.
* **UTW**: Executable generated for the basic implementation.
* **UTWFF**: Executable generated for the implementation with FastFlow.
* **UTWOMP**: Executable generated for the implementation with OpenMP.
* **UTWMPI**: Executable generated for implementation with MPI.