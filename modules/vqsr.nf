//Build SNP Recalibration Model

process SNP_RECALILB {
    input:
        path ref_genome
        path vcf
        path hapmap
        path omni
        path population_g
        path dbsnp
    
    output:
        tuple val sample_id, path("${sample_id}.snf.recal"),path("${sample_id}.snf.tranches"),  emit: recal_model
        path("${sample_id}.snp.plots.R")

    script:
    """
    gatk VariantRecalibrator \
        -R ../reference/${ref_genome} \
        -V ${vcf} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${population_g} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an MQ \
        -mode SNP \
        --output ${sample_id}.snp.recal \
        --tranches-file ${sample_id}.snp.tranches \
        --rscript-file ${sample_id}.snp.plots.R`
    """
}
