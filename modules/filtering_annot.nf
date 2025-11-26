/*
   FIltering & Annotation Nextflow Workflow
========================================================================================
*/

nextflow.enable.dsl=2

// Pipeline Input parameters
params.outdir = '/mnt/d/BI_prj/wgs_proj/hcm'
params.genome = "/mnt/d/BI_prj/wgs_proj/hcm/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.vcf = "/mnt/d/BI_prj/wgs_proj/hcm/HCM.raw_variants.vcf.gz"  // Pattern for FASTQ files
params.dbsnp ="/mnt/d/BI_prj/wgs_proj/hcm/*.dbsnp138.vcf"
params.hapamap = "/mnt/d/BI_prj/wgs_proj/hcm/*.interval_list"
params.omni = "/mnt/d/BI_prj/wgs_proj/hcm/*.interval_list"
params.population_g = "/mnt/d/BI_prj/wgs_proj/hcm/*.interval_list"
/*
   Processes
*/

// SNP filtering
process SNP_FILTER {
    tag "$sample_id"
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
        path ref_genome
        tuple val sample_id, path(vcf)
    
    output:
        tuple val sample_id, 
        path("${sample_id}_snps.vcf.gz"),
        path("${sample_id}.hard_filtered_snps.vcf.gz")

    script:
    """
    gatk SelectVariants \
        -R ${ref_genome} \
        -V ${vcf} \
        -select-type-to-include SNP \
        -O ${sample_id}_snps.vcf.gz
    

    gatk VariantFiltration \
        -R ${ref_genome} \
        -V ${vcf} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -O ${sample_id}.hard_filtered_snps.vcf.gz
    """
}

// INDEL filtering
process INDEL_FILTER {
    tag "$sample_id"
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
        path ref_genome
        tuple val sample_id, path(vcf)
    
    output:
        tuple val sample_id, 
        path("${sample_id}_indel.vcf.gz"),
        path("${sample_id}.hard_filtered_indel.vcf.gz")

    script:
    """
    gatk SelectVariants \
        -V ${vcf} \
        -select-type-to-include INDEL \
        -O ${sample_id}_indel.vcf.gz
    

    gatk VariantFiltration \
        -R ${ref_genome} \
        -V ${vcf} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -filter "SOR > 10.0" --filter-name "SOR10" \
        -O ${sample_id}.hard_filtered_indel.vcf.gz
    """
}

// Assessing filter performance
process ASSESS_FILTER{
    input:
        path vcf
    
    output:
        path "hard_filter_stats.txt"
}

// Extract only high-quality PASS variants
process HIGH_QUAL_FILTER {

    script:
    """
    bcftools view -f PASS cohort\_two\_samples.snpIndelRecal.vcf.gz \
 
    -Oz -o ${sample_id}_high_quality.vcf.gz
 
    bcftools index ${sample_id}_high_quality.vcf.gz
    """
}

workflow {
    // Channel for joint vcf file
    vcf_ch = Channel.fromPath(params.vcf)
                    . map { file ->
                        def sample_id = file.simpleName
                        tuple(sample_id,file)
                    }
    
    ref_ch = Channel.fromPath(params.genome)
    dbsnp_ch = Channel.fromPath(pa)
    
    SNP_FILTER()
    
    
}