#!/bin/bash
# Slurm Script Input Variables
#SBATCH --error=job_out/%x-%A_%a.out
#SBATCH --output=job_out/%x-%A_%a.out
#SBATCH --time=05:00:00
#SBATCH --mem=100GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load lib/fftw/3.3.9-gnu
module load R/4.2.0

# Setup
INPUT_PATH=$(find 02_GeneSetWeights -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Command
start=$SECONDS
Rscript SKAT.R --pathway_data=GeneSets.csv \
				--gene_data=path/to/gencode.v26.GRCh38.genes.gtf \
				--genotype_data=$1 \
				--results_prefix=$2 \
				--alpha_val=0.05 \
				--cisBP=1000000 \
				--phenotype_type=D \
				--weight_type=b \
				--temp_folder=/scratch/${SLURM_JOBID} \
				--weight_directory=$INPUT_PATH \
				--output_dir=03_SKATResults
rm -rf /scratch/${SLURM_JOBID}
end=$SECONDS
echo "duration: $((end-start)) seconds."