# Parallel 3D Volume Analysis

## Overview

This project implements a parallel solution in C using MPI for analyzing time-series data over a 3D spatial grid. The main objective is to compute the following at each time step across a 3D volume:

1. Count of local minima
2. Count of local maxima
3. Global minimum value
4. Global maximum value

This computation is performed in parallel using 3D domain decomposition to enhance performance on large datasets.

## Features

- Parallel processing using MPI
- 3D spatial decomposition based on user-defined PX, PY, PZ process grid
- Efficient communication using non-blocking `MPI_Isend` and `MPI_Irecv`
- Local and global extrema detection for each time step
- Optimized file reading and data distribution
- Timed execution phases (read, compute, total time)

## Project Structure

- `src.c` – Main source code
- `job.sh` – SLURM job script for running the MPI program
- `Group19.pdf` – Project report with implementation details and optimization strategies
- `Assignment.pdf` – Problem statement and test cases

## Compilation

To compile the code:
```bash
mpicc -o executable_file_name src.c
```

## Execution
```bash
sbatch job.sh
```
## Optimization Highlights
- Replaced one-to-many MPI send with MPI_Scatter for better distribution efficiency
- Optimized file I/O by using rank 0 to read and distribute data
- Used non-blocking communication to overlap computation and data exchange
- Avoided unnecessary MPI_Comm splits to reduce overhead
