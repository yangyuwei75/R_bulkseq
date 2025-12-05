#!/bin/bash

INPUT_DIR="Star_output"
OUTPUT_DIR="Gene_Counts"

mkdir -p $OUTPUT_DIR

for tab_file in ${INPUT_DIR}/*_ReadsPerGene.out.tab
do
    filename=$(basename "$tab_file")
    sample_name=${filename%_ReadsPerGene.out.tab}
    tail -n +5 "$tab_file" | cut -f 1,2 > "${OUTPUT_DIR}/${sample_name}.txt"
done
