#!/bin/bash
# Slurm Script Input Variables
#SBATCH --job-name=sim_ref_weights
#SBATCH --error=job_out/%x-%A_%a.out
#SBATCH --output=job_out/%x-%A_%a.out
##SBATCH --partition=cpu2021-bf24
#SBATCH --time=3-00:00:00
#SBATCH --mem=2GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# get the number of files to process per job from the input arguments
FILES_PER_JOB=$1

# calculate the start and end indices using the SLURM_ARRAY_TASK_ID variable and the FILES_PER_JOB input argument
START_INDEX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $FILES_PER_JOB + 1 ))
END_INDEX=$(( $START_INDEX + $FILES_PER_JOB - 1 ))

# loop over the range of indices and process each file
for (( i=$START_INDEX; i<=$END_INDEX; i++ ))
do
    # Setup
    INPUT_PATH=$(find Simulations_Reduced/ -name "*.tsv" | sed -n ${i}p)
    filename=$(basename -- "$INPUT_PATH")
    filename="${filename%.*}"
    pathwayname=$(echo $filename | cut -d'_' -f1)

    subdir="$(dirname $INPUT_PATH)"
    dir1="${subdir##*/}"
    subdir="$(dirname $subdir)"
    dir2="${subdir##*/}"

    # Command
    start=$SECONDS
    java -Djava.io.tmpdir=/scratch/${SLURM_JOBID} -jar ./Jawamix5.jar emmax_select_regions -ig path/to/GTEx.EUR.common.imputed.filtered.hdf5 \
                                                                                                                        -ip ${INPUT_PATH} \
                                                                                                                        -rf SNPAvailability/${pathwayname}.tsv \
                                                                                                                        -o Simulations_Weights/${dir2}/${dir1}/${filename} \
                                                                                                                        -ik path/to/GTEx.EUR.common.imputed.filtered.rescaled.IBS \
                                                                                                                        -p 100000
    rm -rf /scratch/${SLURM_JOBID}
    end=$SECONDS
    echo "duration: $((end-start)) seconds."
done