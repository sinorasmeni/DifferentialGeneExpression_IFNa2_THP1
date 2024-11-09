#!/bin/bash

# Salmon script for quantifying transcript abundances

# Define the directory containing the FASTQ files
INPUT_DIR="path/to/your/fastq/files"
# Define the output directory for Salmon quantification results
OUTPUT_DIR="path/to/output/directory"
# Define the path to the reference transcriptome index
TRANSCRIPTOME_INDEX="path/to/salmon/index"

# Loop through each FASTQ file in the input directory and run Salmon
for file in "$INPUT_DIR"/*.fastq
do
  # Extract the sample name (assuming the filename structure is sample.fastq)
  SAMPLE_NAME=$(basename "$file" .fastq)

  echo "Quantifying $file..."
  
  salmon quant -i "$TRANSCRIPTOME_INDEX" -l A  -r "$file" --validateMappings --numBootstraps 0 --seqBias  -o "$OUTPUT_DIR/$SAMPLE_NAME"_quant
done

