#!/bin/sh
#SBATCH --job-name=fastqc

# set up job environment
module load Java/11.0.27
module load FastQC/0.12.1-Java-11
# working folder: /scratch/brussel/108/vsc10841/RHome

for i in ./Data_monocyte_infiltration/*.fastq.gz
do
        fastqc -o ../QC_RAW/ $i
done
