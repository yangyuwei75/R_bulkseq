#!/bin/bash
#SBATCH --job-name=Mapping_STAR
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# setup job environment
module load STAR

GENOME_DIR="./STAR_Index_hg38_81"
INPUT_DIR="Trimmomatic"
OUTPUT_DIR="Star_output"
CPU=8

mkdir -p $OUTPUT_DIR

for fastq_file in ${INPUT_DIR}/Final_*_R1.fastq.gz
do

    filename=$(basename "$fastq_file")
    temp_name=${filename#Final_}
    sample_name=${temp_name%_R1.fastq.gz}


    STAR --runThreadN $CPU \
         --genomeDir $GENOME_DIR \
         --readFilesIn "$fastq_file" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${OUTPUT_DIR}/${sample_name}_" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --outSAMunmapped Within \
         --outSAMattributes Standard

done
