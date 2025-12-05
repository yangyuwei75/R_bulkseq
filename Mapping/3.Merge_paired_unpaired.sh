#!/bin/sh
#SBATCH --job-name=Trimmomatic

INPUT_DIR="Trimmomatic"
OUTPUT_DIR="Trimmomatic"

for paired_r1 in ${INPUT_DIR}/cutadapt_*_R1.fastq.gz
do
    [ -e "$paired_r1" ] || continue
    
    unpaired_r1="${paired_r1%_R1.fastq.gz}_R1_unpaired.fastq.gz"
    filename=$(basename "$paired_r1")
    temp_name=${filename#cutadapt_}
    sample_name=${temp_name%_R1.fastq.gz}
    final_output="${OUTPUT_DIR}/Final_cutadapt_${sample_name}_R1.fastq.gz"
    cat "$paired_r1" "$unpaired_r1" > "$final_output"

done
