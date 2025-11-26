/*
   Preprocess Module
========================================================================================
*/

// Download SRA file 
process DWNLD_SRA {
    tag "$srr_id"
    publishDir params.outdir, mode: 'copy'
    label 'high_process'

    input:
        val srr_id

    output:
        tuple val(srr_id), path("${srr_id}_1.fastq.gz"), path("${srr_id}_2.fastq.gz"), emit: reads

    script:
    """
    echo "Starting download of $srr_id at \$(date)"
    prefetch $srr_id
    fasterq-dump --split-files --threads 4 $srr_id
    
    # Compress files
    gzip ${srr_id}_1.fastq
    gzip ${srr_id}_2.fastq
    
    echo "Download complete at \$(date)"
    ls -lh ${srr_id}_*.fastq.gz
    """
}

// Quality check FASTQ files
process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/results/fastqc", mode: 'copy'
    label 'low_process'

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        path "*_fastqc.{zip,html}", emit: fastqc_results
        tuple val(sample_id), path(read1), path(read2), emit: reads

    script:
    """
    echo "Running FastQC on ${sample_id}"
    fastqc -t 2 ${read1} ${read2}
    echo "FastQC complete for ${sample_id}"
    """
}

// MultiQC to aggregate results (optional)
process MULTIQC {
    publishDir "${params.outdir}/results/multiqc", mode: 'copy'
    label 'low_process'

    input:
        path fastqc_files

    output:
        path "multiqc_report.html", emit: report
        path "multiqc_data", emit: data

    script:
    """
    multiqc .
    """
}

// Quality triming
process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/results/trimmed_files", mode: 'copy'
    label 'mid_process'

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        tuple val(sample_id), path("${sample_id}_1_trimmed.fastq.gz"), path("${sample_id}_2_trimmed.fastq.gz"), emit: trimmed_reads
        path "${sample_id}_trimmed_report.html", emit: html_report
    
    script:
    """
    echo "Trimming started for ${sample_id}"

    fastp -i ${read1} -o ${sample_id}_1_trimmed.fastq.gz \
          -I ${read2} -O ${sample_id}_2_trimmed.fastq.gz \
          -h ${sample_id}_trimmed_report.html \
          -w 8
    
    echo "Trimming complete for ${sample_id}"
    """
}

// Quality check FASTQ files
process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/results/fastqc", mode: 'copy'
    label 'low_process'

    input:
        tuple val(sample_id), path(read1), path(read2)

    output:
        path "*_fastqc.{zip,html}", emit: fastqc_results
        tuple val(sample_id), path(read1), path(read2), emit: trimmed_reads

    script:
    """
    echo "Running FastQC on trimmed ${sample_id}"
    fastqc -t 2 ${read1} ${read2}
    echo "FastQC complete for ${sample_id}"
    """
}

process MULTIQC_TRIMMED {
    publishDir "${params.outdir}/results/multiqc", mode: 'copy'
    label 'low_process'

    input:
        path fastqc_files

    output:
        path "multiqc_report_trimmed.html", emit: report
        path "multiqc_data_trimmed", emit: data

    script:
    """
    multiqc .
    """
}
