#!/bin/bash

# get the number of files in the input directory
NUM_FILES=$(find 01_DimReduc/ -name "*.tsv" | wc -l)

# submit the SLURM script with the correct number of tasks and arguments
sbatch --array=1-$NUM_FILES ReferenceWeights.sh