#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node=8
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --time=00:20:00         ## wall-clock time limit
#SBATCH --partition=standard    ## can be "standard" or "cpu"

echo `date`
mpirun -np P ./executable_file_name input_data_file_name PX PY PZ NX NY NZ NC output_file_name
echo `date`
