# DifferentialGeneExpression_IFNa2_THP1

## Project Overview 
This project identifies and characterises differentially expressed genes and transcripts in THP-1 cells following a 24-hour stimulation with type I IFN. 

## Objectives
- Quantify transcript and gene expression from RNA-seq data using **Salmon**.
- Identify differentially expressed genes and transcripts with **edgeR**.
- Analyze differential transcript usage with **IsoformSwitchAnalyzeR**.

## Repository Structure
- `quality_control/`: Scripts for running FastQC for data quality checks.
- `trimming/`: Scripts for adapter and quality trimming using TrimGalore.
- `quantification/`: Scripts for expression quantification using Salmon.
- `dge_analysis/`: Scripts for differential gene expression analysis with edgeR.
- `transcript_usage/`: Scripts for transcript usage analysis with IsoformSwitchAnalyzeR.

## Requirements
- **R** (v4.4.1) with packages `edgeR`, `IsoformSwitchAnalyzeR`
- **FastQC**, **TrimGalore**, **Salmon**

## Usage
This analysis pipeline involves the following steps:
1. **Quality Control** with FastQC
2. **Adapter and Quality Trimming** with TrimGalore
3. **Expression Quantification** using Salmon
4. **Differential Expression Analysis** with edgeR
5. **Differential Transcript Usage Analysis** using IsoformSwitchAnalyzeR
