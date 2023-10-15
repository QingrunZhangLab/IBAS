#!/bin/bash
# Slurm Script Input Variables
#SBATCH --error=job_out/%x-%A_%a.out
#SBATCH --output=job_out/%x-%A_%a.out
##SBATCH --partition=cpu2021-bf24
#SBATCH --time=05:00:00
#SBATCH --mem=4GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load R/4.2.0

# Setup
INPUT_PATH=$(find $1 -name "*named.residuals.csv"  | sed -n ${SLURM_ARRAY_TASK_ID}p)
subdir="$(dirname $INPUT_PATH)"
dir1="${subdir##*/}"
subdir="$(dirname $subdir)"
dir2="${subdir##*/}"

# Command
start=$SECONDS
Rscript 01_DimReduc.R --expr_data=${INPUT_PATH} \
				--pathway_data=GeneSets.csv \
				--gene_data=path/to/gencode.v26.GRCh38.genes.gtf \
				--output_dir=./Simulations_Reduced/${dir2}/${dir1}
rm -rf /scratch/${SLURM_JOBID}
end=$SECONDS
echo "duration: $((end-start)) seconds."