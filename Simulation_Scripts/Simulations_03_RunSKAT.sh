#!/bin/bash

# get the number of files in the input directory
NUM_DIRS=$(find Simulations_Weights/ -mindepth 3 -maxdepth 3 -type d | wc -l)

# set the maximum number of jobs to submit at a time
MAX_JOBS=4000

# calculate the number of files to process per job
DIRS_PER_JOB=$(( ($NUM_DIRS + $MAX_JOBS - 1) / $MAX_JOBS ))

# calculate the number of batches
NUM_BATCHES=$(( ($NUM_DIRS + $DIRS_PER_JOB - 1) / $DIRS_PER_JOB ))

echo "Submitting ${NUM_BATCHES} jobs with ${DIRS_PER_JOB} directories to be processed by each job."

# submit the SLURM script as a job array
sbatch --array=1-$NUM_BATCHES Simulation_SKAT.sh $DIRS_PER_JOB $1 $2