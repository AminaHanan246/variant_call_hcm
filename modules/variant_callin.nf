/*
   Variant-Calling Module
========================================================================================
*/

// Index generation
process IDX_GENERATE {
    publishDir "${params.outdir}/results/reference", mode: 'copy'
    label 'mid_process'

    input:
        path ref_genome

    output:
        path "${ref_genome}.fai", emit: ref_idx
        
    script:
    """
    echo "Creating ${ref_genome.baseName} index"

    samtools faidx ${ref_genome}

    echo "Created ${ref_genome.baseName} index"
    """
}


// Dictionary generation
process DICT_GENERATE {
    publishDir "${params.outdir}/results/reference", mode: 'copy'
    label 'low_process'

    input:
        path ref_genome

    output:
        path "${ref_genome.baseName}.dict", emit: ref_dict
        
    script:
    """
    echo "Creating ${ref_genome.baseName} dictionary"

    gatk CreateSequenceDictionary \
        -R ${ref_genome} \
        -O ${ref_genome.baseName}.dict

    echo "Created ${ref_genome.baseName} dictionary"
    """
}


// BWA-MEM indexing
process BWA_IDX {
    publishDir "${params.outdir}/results/reference", mode: 'copy'
    label 'high_process'

    input:
        path ref_genome

    output:
        tuple path(ref_genome), path("${ref_genome}.*"), emit: bwa_index

    script:
    """
    echo "Running BWA index on ${ref_genome.baseName}"

    bwa index ${ref_genome}

    echo "BWA index complete for ${ref_genome.baseName}"

    ls -lh ${ref_genome}.*

    """
}

// BWA-MEM align
process BWA_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'
    label 'high_process'

    input:
        tuple path(ref_genome), path("*")
        tuple val(sample_id), path(read1), path(read2)
        
    output:
        tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: aligned_bam
        path "${sample_id}_flagstat.txt", emit: flagstat
        path "${sample_id}_bwa.log", optional: true

    script:
    """
    echo "Aligning ${sample_id} to reference genome"

    bwa mem -t 8 \
        ${ref_genome} ${read1} ${read2} | samtools sort -@ 8 -o ${sample_id}_sorted.bam -
    
    samtools index ${sample_id}_sorted.bam

    samtools flagstat ${sample_id}_sorted.bam > ${sample_id}_flagstat.txt

    if [ ! -s ${sample_id}_sorted.bam ]; then
    echo "ERROR: BAM file not created"
    exit 1
    fi

    echo "Alignment complete for ${sample_id}"

    ls -lh ${sample_id}_sorted.bam*

    """
}

// Mark duplicate
process MARK_DUPLICATE{
    tag "$sample_id"
    publishDir "${params.outdir}/results/dedup", mode: 'copy'
    label 'mid_process'

    input:
        tuple val(sample_id), path(bam), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}_marked_duplicates.bam"), path("${sample_id}_marked_duplicates.bam.bai"), emit: marked_reads
        path "${sample_id}_duplicate_metrics.txt", emit: metrics  

    script:  
    """
    picard MarkDuplicates \
    -I ${bam} \
    -O ${sample_id}_marked_duplicates.bam \
    -M ${sample_id}_duplicate_metrics.txt \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY LENIENT
    
    echo "=== Duplicate Metrics ==="

    head -8 ${sample_id}_duplicate_metrics.txt | tail -2
    """
}

// Due to better precision of sequencing technologies, BQSR maybe skipped

// process BQSR_TABLE {
//     tag "$sample_id"
//     publishDir "${params.outdir}/results", mode: 'copy'

//     input:
//         tuple val(sample_id), path(bam), path(bai)
//         path ref_genome
//         path ref_fai
//         path ref_dict
//         path dbsnp
//         path known_indels
//         path known_indels2

//     output:
//         tuple val(sample_id), path "${sample_id}_recal_data.table", emit: bqsr_table  

//     script:  
//     """
//     gatk --java-options "-Xmx8G" BaseRecalibrator \
//         --input ${bam} \
//         --reference ${ref_genome} \
//         --known-sites ../resources/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
//         --known-sites ../resources/hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
//         --known-sites ../resources/hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
//         --output ${sample_id}_recal_data.table \
//     """
// }

// process BQSR_APPLY {
//     tag "$sample_id"
//     publishDir "${params.outdir}", mode: 'copy'

//     input:
//         tuple val(sample_id), path(bam), path(bai)
//         path ref_genome
//         path recal_table

//     output:
//         tuple val(sample_id), path("${sample_id}_analysis_ready.bam"), path("${sample_id}_analysis_ready.bam.bai") emit: analysis_reads

//     script:  
//     """
//     gatk --java-options "-Xmx8G" ApplyBQSR \
//         --input ${bam} \
//         --reference ${ref_genome} \
//         --bqsr-recal-file results/${recal_table} \
//         --output ${sample_id}_analysis_ready.bam \
//     """
// }


process HAPLOTYPECALLER {
    tag "$sample_id"
    publishDir "/results/vcf", mode: 'copy'
    label 'high_process'

    input:
        tuple val(sample_id), path(bam), path(bai)
        path ref_genome
        path dbsnp
        path interval_list

    output:
        tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: gVCF
    
    script:  
    """
    gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4" HaplotypeCaller \
    --reference ${ref_genome} \
    --input ${bam} \
    --output ${sample_id}.g.vcf.gz \
    --emit-ref-confidence GVCF \
    --dbsnp ${dbsnp} \
    --intervals ${interval_list} \
    --native-pair-hmm-threads 4 \
    --standard-min-confidence-threshold-for-calling 10.0 \
    --annotation QualByDepth \
    --annotation Coverage \
    --annotation FisherStrand \
    --annotation StrandOddsRatio \
    --annotation MappingQualityRankSumTest \
    --annotation ReadPosRankSumTest \
    --annotation RMSMappingQuality
    """
}

process COMBINE_GVCF {
    publishDir "${params.outdir}/results/vcf", mode: 'copy'
    label 'mid_process'

    input:
        path vcf_list
        path ref_genome

    output:
        tuple path("HCM_combined.g.vcf.gz"), path("HCM_combined.g.vcf.gz.tbi"), emit:combine_vcf

    script:
    """
    gatk CombineGVCFs \
        -R ${ref_genome} \
        ${vcf_list.collect{ "-V ${it}" }.join(" ")} \
        -O HCM_combined.g.vcf.gz
    """
}


process JOINT_GENOTYPING {
    publishDir "${params.outdir}/results/vcf", mode: 'copy'
    label 'id_process'

    input:
        tuple path(vcf), path(vcf_idx)
        path ref_genome

    output:
        tuple path("HCM.raw_variants.vcf.gz"), path("HCM.raw_variants.vcf.gz.tbi"), emit:raw_vcf

    script:
    """
     gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R ${ref_genome} \
        -V ${vcf} \
        -O HCM.raw_variants.vcf.gz
    """
}

process ANALYSE_DATA {
    publishDir "${params.outdir}/results/vcf", mode: 'copy'
    label 'low_process'

    input:
        tuple path(vcf), path(vcf_idx)

    output:
        path "variant_metrics.table", emit:variant_metrics

    script:
    """
    gatk VariantsToTable \
        -V ${vcf} \
        -F CHROM -F POS -F TYPE -F QD -F FS -F MQ -F SOR -F MQRankSum -F ReadPosRankSum \
        -O variant_metrics.table

    """
}