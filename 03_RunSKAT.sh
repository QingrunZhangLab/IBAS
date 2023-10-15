#!/bin/bash

# get the number of files in the input directory
NUM_FILES=$(find 02_GeneSetWeights -maxdepth 1 -type d | wc -l)

# submit the SLURM script with the correct number of tasks and arguments
sbatch --array=2-$NUM_FILES --job-name=skat_$2 SKAT.sh $1 $2