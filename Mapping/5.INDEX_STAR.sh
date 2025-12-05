#!/bin/bash
#SBATCH --mem=50G      
#SBATCH --cpus-per-task=8
module load STAR

GTF_FILE="./Reference_genome/Homo_sapiens.GRCh38.81.gtf" 
FASTA_FILE="./Reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENOME_DIR="./STAR_Index_hg38_81"

OVERHANG=100

STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir $GENOME_DIR \
    --genomeFastaFiles $FASTA_FILE \
    --sjdbGTFfile $GTF_FILE \
    --sjdbOverhang $OVERHANG
