#!/bin/bash
# Slurm Script Input Variables
#SBATCH --job-name=sim_skat_runs
#SBATCH --error=job_out/%x-%A_%a.out
#SBATCH --output=job_out/%x-%A_%a.out
##SBATCH --partition=cpu2021-bf24
#SBATCH --time=3-00:00:00
#SBATCH --mem=100GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load lib/fftw/3.3.9-gnu
module load R/4.2.0

# get the number of files to process per job from the input arguments
DIRS_PER_JOB=$1

# calculate the start and end indices using the SLURM_ARRAY_TASK_ID variable and the DIRS_PER_JOB input argument
START_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $DIRS_PER_JOB + 1 ))
END_INDEX=$(( $START_INDEX + $DIRS_PER_JOB - 1 ))

# loop over the range of indices and process each file
for (( i=$START_INDEX; i<=$END_INDEX; i++ ))
do
    # Setup
    INPUT_PATH=$(find Simulations_Weights/ -mindepth 3 -maxdepth 3 -type d | sed -n ${i}p)
    filename=$(basename -- "$INPUT_PATH")
    filename="${filename%.*}"
    pathwayname=$(echo $filename | cut -d'_' -f1)

    subdir="$(dirname $INPUT_PATH)"
    dir1="${subdir##*/}"
    subdir="$(dirname $subdir)"
    dir2="${subdir##*/}"

    # Command
    start=$SECONDS
	Rscript SKAT.R --pathway_data=GeneSets.csv \
					--gene_data=path/to/gencode.v26.GRCh38.genes.gtf \
					--genotype_data=path/to/split_plink/T1D \
					--results_prefix=$2 \
					--alpha_val=0.05 \
					--cisBP=1000000 \
					--phenotype_type=D \
					--weight_type=$3 \
					--temp_folder=/scratch/${SLURM_JOBID} \
					--weight_directory=$INPUT_PATH \
					--output_dir=Simulations_SKATResults/${dir2}/${dir1}
	rm -rf /scratch/${SLURM_JOBID}
    end=$SECONDS
    echo "duration: $((end-start)) seconds."
done