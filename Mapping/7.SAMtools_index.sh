#!/bin/bash
#SBATCH --cpus-per-task=4

module load SAMtools/1.18-GCC-12.3.0
INPUT_DIR="Star_output"

for bam_file in ${INPUT_DIR}/*Aligned.sortedByCoord.out.bam
do
    samtools index -@ 4 "$bam_file"
done
