/*
   Variant-Calling Nextflow Workflow
========================================================================================
*/
include { DWNLD_SRA; FASTQC; FASTP; FASTQC_TRIMMED } from './modules/preprocess.nf'
include { IDX_GENERATE; DICT_GENERATE; BWA_IDX; BWA_ALIGN; MARK_DUPLICATE; HAPLOTYPECALLER; COMBINE_GVCF; JOINT_GENOTYPING; ANALYSE_DATA } from './modules/variant_callin.nf'

workflow {
    // Channel for reference genome
    ref_ch = Channel.fromPath(params.genome)
    
    // Downloaded FASTQ files
    if (params.reads) {
        reads_ch = Channel
            .fromFilePairs(params.reads, size: 2)
            .map { sample_id, files -> 
                tuple(sample_id, files[0], files[1])
            }
    }
    
    // SRA files download
    else if (params.srr) {
        if (file(params.srr).exists()) {
            srr_ch = Channel
                .fromPath(params.srr)
                .splitText()
                .map { it.trim() }
        } 
        else {
            srr_ch = Channel.of(params.srr) 
        }
        
        // Download SRA files
        DWNLD_SRA(srr_ch)
        reads_ch = DWNLD_SRA.out.reads
    }

    else {
        error "Please provide either --reads (for existing FASTQ files) or --srr (to download from SRA)"
    }
        
    // Run FASTQC on existing files
    FASTQC(reads_ch)
        
    // Trimming reads
    FASTP(FASTQC.out.reads)

    // Run FASTQC on trimmed reads
    FASTQC_TRIMMED(FASTP.out.trimmed_reads)

    // // Run MultiQC to aggregate trimmed results
    // MULTIQC_TRIMMED(FASTQC_TRIMMED.out.fastqc_results.collect())

    // Run samtools FASTA indexing
    IDX_GENERATE(ref_ch)

    // Create sequence dictionart
    DICT_GENERATE(ref_ch)

    // BWA index 
    BWA_IDX(ref_ch)

    // BWA Align
    BWA_ALIGN(BWA_IDX.out.bwa_index, FASTP.out.trimmed_reads)

    // // Group all QC files
    // qc_files = FASTQC.out.fastqc_results
    //     .mix(FASTP.out.html_report)
    //     .mix(FASTQC_TRIMMED.out.fastqc_results)
    //     .mix(BWA_ALIGN.out.flagstat)
    //     .collect()
    
    // MULTIQC(qc_files)

    // Mark duplicates
    MARK_DUPLICATE(BWA_ALIGN.out.aligned_bam)

    dbsnp_ch = Channel.fromPath(params.dbsnp)
    interval_ch = Channel.fromPath(params.interval_list)

    // Haplotype Caller
    HAPLOTYPECALLER(
        MARK_DUPLICATE.out.marked_reads,
        ref_ch,
        dbsnp_ch,
        interval_ch 
        )
    
    // Creating list of all VCF files generated
    HAPLOTYPECALLER.out.gVCF
            .map { sample_id, vcf, tbi -> vcf }
            .collect()
            .set{vcf_list}
    
    
    COMBINE_GVCF(vcf_list,ref_ch)

    JOINT_GENOTYPING(COMBINE_GVCF.out.combine_vcf,ref_ch)

    ANALYSE_DATA(JOINT_GENOTYPING.out.raw_vcf)
    
}