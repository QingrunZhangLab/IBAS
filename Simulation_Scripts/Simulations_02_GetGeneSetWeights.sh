#!/bin/bash

# get the number of files in the input directory
NUM_FILES=$(find Simulations_Reduced/ -name "*.tsv" | wc -l)

# set the maximum number of jobs to submit at a time
MAX_JOBS=4000

# calculate the number of files to process per job
FILES_PER_JOB=$(( ($NUM_FILES + $MAX_JOBS - 1) / $MAX_JOBS ))

# calculate the number of batches
NUM_BATCHES=$(( ($NUM_FILES + $FILES_PER_JOB - 1) / $FILES_PER_JOB ))

echo "Submitting ${NUM_BATCHES} jobs with ${FILES_PER_JOB} files to be processed by each job."

# submit the SLURM script as a job array
sbatch --array=1-$NUM_BATCHES Simulation_ReferenceWeights.sh $FILES_PER_JOB