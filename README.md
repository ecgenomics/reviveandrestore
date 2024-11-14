# Pipeline Overview

This Nextflow script performs a DNA sequencing data processing pipeline tailored for capturing and isolating endogenous or specific DNA (hDNA) regions in human samples. The pipeline takes raw sequencing reads and processes them through several steps to produce high-quality, target-enriched BAM files ready for downstream analyses. Here’s a breakdown of the main steps:
	1.	Demultiplexing: Raw reads from pooled sequencing data are separated by sample-specific barcodes, enabling differentiation of multiple samples within a single sequencing run.
	2.	Trimming: Adapter sequences and low-quality bases are removed from the reads, resulting in cleaner data. This step generates trimmed paired-end files that are essential for accurate mapping.
	3.	Alignment: Trimmed reads are aligned to the human reference genome (hg19) using BWA. Both single-end and paired-end alignments are conducted, depending on the sequencing strategy, producing BAM files for each.
	4.	Merging & Deduplication: BAM files from single- and paired-end reads are merged, and duplicate reads are removed. Deduplication is crucial for reducing PCR artifacts and obtaining a more accurate representation of the DNA sample.
	5.	Quality Filtering: Aligned reads are filtered by quality score, retaining only high-confidence reads. This step enhances the reliability of downstream analyses by excluding low-quality data.
	6.	On-Target Filtering: Using a BED file specifying genomic regions of interest, only reads that overlap with these target regions are retained. This on-target enrichment focuses the data on relevant areas, increasing the signal for downstream capture and variant analysis.

The final output consists of high-quality, deduplicated, and on-target BAM files that are optimized for studies requiring targeted DNA capture, such as studies on genetic variation in specific genomic regions. This pipeline automates several steps essential for obtaining reliable, high-quality data for capture-based sequencing applications.

# Requirements

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
