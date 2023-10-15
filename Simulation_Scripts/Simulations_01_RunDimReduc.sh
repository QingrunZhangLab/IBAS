#!/bin/bash

peerOutput_dir=path/to/peerOutputs

# get the number of files in the input directory
NUM_FILES=$(find $peerOutput_dir -name "*named.residuals.csv" | wc -l)

# submit the SLURM script with the correct number of tasks and arguments
sbatch --array=1-$NUM_FILES --job-name=sim_dim_reduc Simulation_DimReduc.sh $peerOutput_dir