# Pipeline Overview

This pipeline processes high-throughput sequencing data to quantify endogenous or human DNA, preparing equi-endogenous pools for capture. Based on Mapping Capture Experiments, this workflow performs tasks like demultiplexing, trimming, mapping, and quality filtering to create high-quality data for downstream analyses.

#Requirements

This pipeline assumes access to the following tools:

	•	sabre
	•	picard
	•	bwa
	•	samtools
	•	bedtools

Please ensure paths are correctly set for each tool as specified in the parameters.

#Inputs

	•	params.reads: Specifies the path to input FASTQ files for sequencing reads.
	•	params.genome: Points to the reference genome (e.g., hg19) used for alignment.
	•	params.bedfile: Specifies the BED file for regions of interest, used in target filtering.

# Pipeline Steps

## 0. Multiplexing

The multiplexing process separates samples based on their barcodes using sabre, creating demultiplexed FASTQ files for each sample.

	•	Output: Demultiplexed FASTQ files saved to params.dm.

## 1. Trimming

The trim process removes adapters and performs quality trimming using fastp. This process outputs both paired and unpaired trimmed FASTQ files.

	•	Output: Trimmed FASTQ files and a report in HTML format saved to params.tr_out.

## 2. Alignment

Alignment is divided into two parts:

	•	2.1 Single-End Alignment (align_se): Aligns merged reads (single-end) to the genome using bwa. Output BAM files are indexed and contain read group information.
	•	2.2 Paired-End Alignment (align_pe): Aligns paired-end reads, filtering based on insert size to ensure high-quality mappings.

Both single-end and paired-end BAM files are saved to params.bam_hg19.

## 3. Merging and Duplicate Removal

	•	3.1 Merging (merge_sepe): Merges single- and paired-end BAM files using picard, sorts, and indexes the resulting BAM files.
	•	3.2 Duplicate Removal (rmdup): Removes duplicate reads to avoid bias in quantification, using picard.
	•	Output: Deduplicated BAM files saved to params.bam_hg19_no_dups.

##4. Quality Filtering

The filterQuality process sorts BAM files and filters them based on a minimum mapping quality threshold of 30.

	•	Output: Filtered BAM files are saved to params.BAM_hg19_filtered.

## 5. Target Region Filtering

The ontarget process extracts reads mapping to target regions specified in the BED file, isolating only on-target reads.

	•	Output: On-target BAM files saved to params.BAM_ontarget.

# Workflow Execution

The workflow links each step sequentially, starting from demultiplexing and ending with filtering on-target regions. This ensures that each output is correctly passed to the next step in the pipeline for efficient and reproducible processing of sequencing data.
