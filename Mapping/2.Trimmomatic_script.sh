#!/bin/sh
#SBATCH --job-name=Trimmomatic

# setup job environment
module load Java/17.0.15
module load Trimmomatic
#  working folder: /scratch/brussel/108/vsc10841/RHome

INPUT_DIR="./Data_monocyte_infiltration"
OUTPUT_DIR="Trimmomatic"

for r1_file in ${INPUT_DIR}/*_R1_001.fastq.gz
    do
     
     filename=$(basename "$r1_file")
     base_name=${filename%_R1_001.fastq.gz}
     r2_file="${INPUT_DIR}/${base_name}_R2_001.fastq.gz"
    
     java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.40.jar PE \
     "$r1_file" "$r2_file" \
     "${OUTPUT_DIR}/cutadapt_${base_name}_R1.fastq.gz" \
     "${OUTPUT_DIR}/cutadapt_${base_name}_R1_unpaired.fastq.gz" \
     "${OUTPUT_DIR}/cutadapt_${base_name}_R2_ignore.fastq.gz" \
     "${OUTPUT_DIR}/cutadapt_${base_name}_R2_unpaired.fastq.gz" \
     ILLUMINACLIP:/scratch/brussel/108/vsc10841/RHome/adpaters_plus_polyA.fa:2:30:10 \
     LEADING:3 \
     TRAILING:3 \
     SLIDINGWINDOW:4:15 \
     MINLEN:36 \
     # ILLUMINACLIP: cut adapter or illumina-specific sequences
     # LEADING, TRAILING: cut bases off the start of a read
     # MINLEN: Drop the read if it is below a specified length.
     # SLIDINGWINDOW: number of base = 4:Quality lower than 15, cut them off
    done
