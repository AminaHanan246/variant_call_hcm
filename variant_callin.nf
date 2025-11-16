/*
   Variant-Calling Nextflow Workflow
========================================================================================
*/

nextflow.enable.dsl=2

// Pipeline Input parameters
params.outdir = '/mnt/d/BI_prj/wgs_proj/hcm'
params.genome = "/mnt/d/BI_prj/wgs_proj/hcm/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.reads = "/mnt/d/BI_prj/wgs_proj/hcm/*_{1,2}.fastq.gz"  // Pattern for FASTQ files
params.srr = null  // Optional: for download mode

/*
   Processes
*/

// Download SRA file (optional - commented out by default)
process DWNLD_SRA {
    tag "$srr_id"
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir params.outdir, mode: 'copy'

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
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}/results/fastqc", mode: 'copy'

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
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}/results/multiqc", mode: 'copy'

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
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir params.outdir, mode: 'copy'

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
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}/results/fastqc", mode: 'copy'

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
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}/results/multiqc", mode: 'copy'

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


// Index generation
process IDX_GENERATE {
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}", mode: 'copy'

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
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}", mode: 'copy'

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
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}", mode: 'copy'

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
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(sample_id), path(read1), path(read2)
        tuple path(ref_genome), path(genome_idx)
        
    output:
        tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: aligned_bam
        path "${sample_id}_flagstat.txt", emit: flagstat

    script:
    """
    echo "Aligning ${sample_id} to reference genome"

    bwa mem -t 8 \\
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
        ${ref_genome} ${read1} ${read2} | \\
    samtools sort -@ 8 -o ${sample_id}_sorted.bam -
    
    samtools index ${sample_id}_sorted.bam
    samtools flagstat ${sample_id}_sorted.bam > ${sample_id}_flagstat.txt

    echo "Alignment complete for ${sample_id}"

    ls -lh ${sample_id}_sorted.bam*

    """
}


workflow {
    // Channel for reference genome
    ref_ch = Channel.fromPath(params.genome)
    
    // Option 1: Use pre-downloaded FASTQ files
    if (params.reads) {
        reads_ch = Channel
            .fromFilePairs(params.reads, size: 2)
            .map { sample_id, files -> 
                tuple(sample_id, files[0], files[1])
            }
    }
    
    // Option 2: Download files first, then process
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
    BWA_ALIGN(FASTP.out.trimmed_reads,BWA_IDX.out.bwa_index)

    // Group all QC files
    qc_files = FASTQC.out.fastqc_results
        .mix(FASTP.out.html_report)
        .mix(FASTQC_TRIMMED.out.fastqc_results)
        .mix(BWA_ALIGN.out.flagstat)
        .collect()
    
    MULTIQC(qc_files)
    
}