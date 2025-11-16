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
        tuple val(sample_id), path("${sample_id}_1_trimmed.fastq"), path("${sample_id}_2_trimmed.fastq"), emit: trimmed_reads
    
    script:
    """
    echo "Trimming started for ${sample_id}"

    fastp -i ${read1} -o ${sample_id}_1_trimmed.fastq \
          -I ${read2} -O ${sample_id}_2_trimmed.fastq \
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

    script:
    """
    echo "Running FastQC on trimmed ${sample_id}"
    fastqc -t 2 ${read1} ${read2}
    echo "FastQC complete for ${sample_id}"
    """
}

// Index generation
process IDX_GENERATE {
    conda "/mnt/d/BI_prj/wgs_proj/hcm/work/conda/"
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path ref_genome

    output:
        path "${ref_genome.baseName}.fai", emit: ref_idx
        
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
        -R Homo_sapiens_assembly38.fasta \
        -O Homo_sapiens_assembly38.dict

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
        path "${ref_genome}.ann", emit: ref_ann
        path "${ref_genome}.bwt", emit: ref_bwt
        path "${ref_genome}.pac", emit: ref_pac
        path "${ref_genome}.sa", emit: ref_sa
        path "${ref_genome}.amb", emit: ref_amb

    script:
    """
    echo "Running BWA index on ${ref_genome}"

    bwa index ${ref_genome}

    echo "BWA index complete for ${ref_genome}"
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
        
    // Run FASTQC on existing files
        FASTQC(reads_ch)
        
    // Run MultiQC to aggregate results
    MULTIQC(FASTQC.out.fastqc_results.collect())

    // Trimming reads
    FASTP(reads_ch)

    // Run FASTQC on trimmed reads
    FASTQC_TRIMMED(FASTP.out.trimmed_reads)

    // Run MultiQC to aggregate trimmed results
    MULTIQC(FASTQC_TRIMMED.out.fastqc_results.collect())

    // Run samtools FASTA indexing
    IDX_GENERATE(ref_ch)

    // Create sequence dictionart
    DICT_GENERATE(ref_ch)

    // BWA index 
    BWA_IDX(ref_ch)

    
    }
    
    else {
        error "Please provide either --reads (for existing FASTQ files) or --srr (to download from SRA)"
    }
}