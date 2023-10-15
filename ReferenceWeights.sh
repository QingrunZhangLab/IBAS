#!/bin/bash
# Slurm Script Input Variables
#SBATCH --job-name=ref_weights
#SBATCH --error=job_out/%x-%A_%a.out
#SBATCH --output=job_out/%x-%A_%a.out
#SBATCH --time=5:00:00
#SBATCH --mem=2GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Setup
INPUT_PATH=$(find 01_DimReduc/ -name "*.tsv" | sed -n ${SLURM_ARRAY_TASK_ID}p)
filename=$(basename -- "$INPUT_PATH")
filename="${filename%.*}"
pathwayname=$(echo $filename | cut -d'_' -f1)

# Command
start=$SECONDS
java -Djava.io.tmpdir=/scratch/${SLURM_JOBID} -jar /path/to/Jawamix5.jar emmax_select_regions -ig /path/to/GTEx.EUR.common.imputed.filtered.hdf5 \
																													-ip ${INPUT_PATH} \
																													-rf SNPAvailability/${pathwayname}.tsv \
																													-o 02_GeneSetWeights/${filename} \
																													-ik /path/to/GTEx.EUR.common.imputed.filtered.rescaled.IBS \
																													-p 100000
rm -rf /scratch/${SLURM_JOBID}
end=$SECONDS
echo "duration: $((end-start)) seconds."