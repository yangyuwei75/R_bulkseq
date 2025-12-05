#!/bin/sh
#SBATCH --job-name=fastqc

# setup job environment
module load Java/11.0.27
module load FastQC/0.12.1-Java-11
# working folder: /scratch/brussel/108/vsc10841/RHome

for i in ./Trimmomatic/Final_cutadapt_*
do
    fastqc -o /scratch/brussel/108/vsc10841/RHome/QC_trimmed/ $i
done
