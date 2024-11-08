#!/bin/bash

#FastQC script for quality control of FASTQ files

# Define the directory containing the FASTQ files
INPUT_DIR="path/to/your/fastq/files"   # Replace with your directory path
# Define the directory to save FastQC results
OUTPUT_DIR="path/to/output/directory"  # Replace with your output path

# Loop through each FASTQ file in the input directory and run FastQC
for file in "$INPUT_DIR"/*.fastq
do
  echo "Running FastQC on $file..."
  fastqc "$file" -o "$OUTPUT_DIR"
done
