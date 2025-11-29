/*
   FIltering & Annotation Nextflow Workflow
========================================================================================
*/

// SNP filtering
process SNP_FILTER {
    
    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
        path ref_genome
        path(vcf)
    
    output:
        path("cohort.snps.vcf.gz"), emit: raw_snps
        path("cohort.hard_filtered_snps.vcf.gz"), emit: filtered_snps

    script:
    """
    gatk SelectVariants \
        -R ${ref_genome} \
        -V ${vcf} \
        -select-type-to-include SNP \
        -O cohort.snps.vcf.gz
    

    gatk VariantFiltration \
        -R ${ref_genome} \
        -V cohort.snps.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -O cohort.hard_filtered_snps.vcf.gz
    """
}

// INDEL filtering
process INDEL_FILTER {

    publishDir "${params.outdir}/filtered_variants", mode: 'copy'

    input:
        path ref_genome
        path(vcf)
    
    output: 
        path("cohort.indel.vcf.gz"), emit: raw_indel
        path("cohort.hard_filtered_indel.vcf.gz"), emit: filtered_indel

    script:
    """
    gatk SelectVariants \
        -R ${ref_genome} \
        -V ${vcf} \
        -select-type-to-include INDEL \
        -O cohort.indel.vcf.gz
    

    gatk VariantFiltration \
        -R ${ref_genome} \
        -V cohort.indel.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -filter "SOR > 10.0" --filter-name "SOR10" \
        -O cohort.hard_filtered_indel.vcf.gz
    """
}

// Assessing filter performance
process ASSESS_FILTER{
    input:
        path vcf
    
    output:
        path "hard_filter_stats.txt", emit: filter_stat
}

process MERGE_VCF {
    input:
        path snp_filtered_vcf
        path indel_filtered_vcf
    output:
        path ("cohort_filtered.vcf.gz"), emit: filtered_vcf
    
    script:
    """
        gatk MergeVcfs \
            -I ${snp_filtered_vcf} \
            -I ${indel_filtered_vcf} \
            -O cohort_filtered.vcf.gz 
    """
}

// Extract only high-quality PASS variants
process HIGH_QUAL_FILTER {
    input:
        path merged_vcf
    
    output:
        path("cohort_high_quality.vcf.gz"), emit: filtered_vcf

    script:
    """
    bcftools view -f PASS ${merged_vcf} \
 
    -Oz -o cohort_high_quality.vcf.gz
 
    bcftools index cohort_high_quality.vcf.gz
    """
}

// 
