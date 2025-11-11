#!/usr/bin/env nextflow

/*
 * Download fasta files
 */
process DWNLD_SRA {
    tag "$srr_id"
    conda "preprocess"

    publishDir params.outdir, mode: 'symlink'

    input:
        val srr_id

    output:
        tuple path("${srr_id}_1.fastq.gz"), path("${srr_id}_2.fastq.gz"), emit: reads

    script:
    """
    fasterq-dump $srr_id --split-files --gzip --threads 4
    """

}
process FASTQC {

    conda "preprocess"
    publishDir "results/", mode: 'symlink'

    input:
    path reads

    output:
    path "${reads.simpleName}_fastqc.zip", emit: zip
    path "${reads.simpleName}_fastqc.html", emit: html

    script:
    """
    fastqc $reads
    """
}

process MULTIQC {

    conda "preprocess"
    publishDir "results/", mode: 'symlink'

    input:
    path 'results/*'
    val output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}

process TRIMMOMATIC {

    conda "preprocess"
    publishDir "*", mode: 'symlink'

    input:
    tuple path(r1), path(r2)

    output:
    path "${r1.simpleName}_trimmed.fastq.gz", emit: trimmed_r1
    path "${r2.simpleName}_trimmed.fastq.gz", emit: trimmed_r2
    path "${r1.simpleName}_trimming_report.txt", emit: trimming_report_r1
    path "${r1.simpleName}_trimmed_fastqc.zip", emit: fastqc_zip_r1
    path "${r1.simpleName}_trimmed_fastqc.html", emit: fastqc_html_r1
    path "${r2.simpleName}_trimming_report.txt", emit: trimming_report_r2
    path "${r2.simpleName}_trimmed_fastqc.zip", emit: fastqc_zip_r2
    path "${r2.simpleName}_trimmed_fastqc.html", emit: fastqc_html_r2

    script:
    """
    trimmomatic PE -threads 4 \\
     $r1 $r2 \\
     ${r1.simpleName}_trimmed.fastq.gz ${r2.simpleName}_trimmed.fastq.gz \\ 
     TRAILING:10 -phred33

    fastqc ${r1.simpleName}_trimmed.fastq.gz
    fastqc ${r2.simpleName}_trimmed.fastq.gz

    echo "Trimming complete {$r1.simpleName}"
    """
}
    