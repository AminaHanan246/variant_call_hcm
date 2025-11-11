#!/usr/bin/env nextflow

// Module INCLUDE statements
include { DWNLD_SRA; FASTQC; MULTIQC; TRIMMOMATIC } from './modules/preprocess.nf'

/*
 * Pipeline parameters
 */

// Primary input
params.srr_ids = "data/srr_id.csv"
params.report_id = "var_call"
params.outdir_downloads = ""


workflow {

    // Create input channel from a file path
    srr_ch = channel
       .fromPath(params.srr_ids)
       .splitText()
	
    // Download Fasta
    DWNLD_SRA(srr_ch)
	
    // Call processes
    FASTQC(DWNLD_SRA.out.reads)
    
    // Aggregate report
    MULTIQC(
        FASTQC.out.zip.mix(
        FASTQC.out.html,
        ).collect(),
        params.report_id
        )
   
    // Trimming PE
    TRIMMOMATIC(DWNLD_SRA.out.reads)
 }