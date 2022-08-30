#!/bin/bash
# Slurm Script Input Variables
#SBATCH --job-name=ref_weights
#SBATCH --array=1-120
#SBATCH --error=job_out/%x-%A_%a.out
#SBATCH --output=job_out/%x-%A_%a.out
#SBATCH --time=5:00:00
#SBATCH --mem=2GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Setup
INPUT_PATH=$(find ReducedPathways/ -name "*.tsv" | sed -n ${SLURM_ARRAY_TASK_ID}p)
filename=$(basename -- "$INPUT_PATH")
filename="${filename%.*}"
pathwayname="${filename#*_}"

# Command
start=$SECONDS
java -jar Jawamix5.jar emmax_select_regions -ig path/to/reference/data/genotypes.hdf5 \
																				-ip ${INPUT_PATH} \
																				-rf SNPAvailability/${pathwayname}.tsv \
																				-o PathwayWeights/${filename} \
																				-ik path/to/reference/data/kinship.rescaled.IBS \
																				-p 100000
end=$SECONDS
echo "duration: $((end-start)) seconds."