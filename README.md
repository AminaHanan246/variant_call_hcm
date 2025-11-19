\# Germline Variant calling using Nextflow

A Nextflow workflow for germline variant calling following GATK best practices, from raw FASTQ files through joint genotyping.



\## Overview

This pipeline processes paired-end sequencing data (WGS/WES) to identify germline variants. It handles everything from quality control through joint genotyping of multiple samples, producing analysis-ready VCF files.



\## Pipeline Steps

1. Quality Control - FastQC on raw reads

2\. Trimming - Fastp for adapter removal and quality trimming

3\. Post-trim QC - FastQC on trimmed reads

4\. Reference Preparation - Index generation and dictionary creation

5\. Alignment - BWA-MEM2 alignment to reference genome

6\. Duplicate Marking - Picard MarkDuplicates

7\. Variant Calling - GATK HaplotypeCaller in GVCF mode

8\. Joint Genotyping - CombineGVCFs and GenotypeGVCFs for multi-sample calling

10\. Variant Metrics - Extract key variant statistics for QC

