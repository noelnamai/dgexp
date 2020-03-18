#!/usr/bin/env nextflow

params.help = null

println """\

D I F F E R E N T I A L  G E N E  E X P R E S I O N  A N A L Y S I S
====================================================================

Samtools   : 1.9
Gffread    : 0.11.8
HISAT2     : 2.1.0
HTSeq      : 0.11.1
DESeq2     : 1.26.0
Start time : $workflow.start
"""
.stripIndent()

/* 
 * parse the input parameters.
 */
gff3   = file(params.gff3)
genome = file(params.genome)

/*
 * read all the fastq.gz files pairs into a channel of the form:
 * [[state1_rep1, [state1_rep1.fastq.gz, state1_rep1_2.fastq.gz]], [state1_rep2, [state1_rep2.fastq.gz, state1_rep2_2.fastq.gz]], ...]
 */
Channel.fromPath(params.reads)
    .ifEmpty{error "Cannot find any reads matching: ${params.reads}"}
    .map{file -> [file.name.substring(6, 17), file]}
    .groupTuple(by: 0, sort: true)
    .set{read_pairs_ch}

/*
 * trim the fastq reads using trimmomatic
 */
process trim_reads {

    cpus = 2
    tag "$state_replicate" 
    container "noelnamai/asimov:1.0"

    input:
    set state_replicate, file(reads) from read_pairs_ch

    output:
    set state_replicate, file("${reads[0].simpleName}.trimmed.fastq.gz"), file("${reads[1].simpleName}.trimmed.fastq.gz") into trimmed_read_pairs_1_ch
    set state_replicate, file("${reads[0].simpleName}.trimmed.fastq.gz"), file("${reads[1].simpleName}.trimmed.fastq.gz") into trimmed_read_pairs_2_ch

    script:
    """
    java -jar /usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
        -phred33 \
        -threads ${task.cpus} \
        ${reads[0]} ${reads[0].simpleName}.trimmed.fastq.gz \
        ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    java -jar /usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
        -phred33 \
        -threads ${task.cpus} \
        ${reads[1]} ${reads[1].simpleName}.trimmed.fastq.gz \
        ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process check_read_quality {

    cpus = 2
    tag "$state_replicate" 
    container "noelnamai/asimov:1.0"

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    set state_replicate, file(read_1), file(read_2) from trimmed_read_pairs_1_ch

    output:
    file("*.html") into fastqc_ch

    script:
    """
    fastqc ${read_1} --outdir . --threads ${task.cpus}
    fastqc ${read_2} --outdir . --threads ${task.cpus}
    """
}

/*
 * convert the annotated gff3 to gtf using gffread
 */
 process convert_gff3_to_gtf {

    cpus = 2
    container "noelnamai/asimov:1.0"

    input:
    file gff3
      
    output:
    file("${gff3.baseName}.gtf") into gtf_1_ch
    file("${gff3.baseName}.gtf") into gtf_2_ch

    script:
    """
    gffread -T ${gff3} -o ${gff3.baseName}.gtf
    """
}

/*
 * extract exons and splice junctions from the gff3 file
 */
process extract_exons_and_ss {

    cpus = 2
    container "noelnamai/asimov:1.0"

    input:
    file(gtf) from gtf_1_ch
      
    output:
    set file("${gtf.baseName}.exons.tsv"), file("${gtf.baseName}.splicesites.tsv") into extracted_exons_and_ss_ch

    script:
    """    
    hisat2_extract_exons.py ${gtf.baseName}.gtf > ${gtf.baseName}.exons.tsv
    hisat2_extract_splice_sites.py ${gtf.baseName}.gtf > ${gtf.baseName}.splicesites.tsv
    """
}

/*
 * build the genome index for mapping using hisat2
 */
process build_genome_index {
    
    cpus = 2
    container "noelnamai/asimov:1.0"

    input:
    file genome
    set file(extracted_exons), file(extracted_ss) from extracted_exons_and_ss_ch
      
    output:
    file("${genome.simpleName}.index.*") into genome_index_ch

    script:
    """
    hisat2-build --exon ${extracted_exons} --ss ${extracted_ss} -p ${task.cpus} ${genome} ${genome.simpleName}.index
    """
}

/*
 * align reads to the reference genome using hisat2
 */
process map_reads_to_reference {
    
    cpus = 2
    tag "$state_replicate"
    container "noelnamai/asimov:1.0"

    input:
    file genome
    file index from genome_index_ch
    set state_replicate, file(read_1), file(read_2) from trimmed_read_pairs_2_ch

    output:
    set state_replicate, file("${state_replicate}.sam") into aligned_sam_ch

    script:
    """
    hisat2 --dta --phred33 -p ${task.cpus} -x ${genome.simpleName}.index -U ${read_1},${read_2} -S ${state_replicate}.sam
    """
}

/*
 * convert sam to sorted bam files using samtools
 */
process convert_sam_to_bam {

    cpus = 2
    tag "$state_replicate"
    container "noelnamai/asimov:1.0"

    input:
    set state_replicate, file(sam_file) from aligned_sam_ch

    output:
    set state_replicate, file("${sam_file.baseName}.bam") into aligned_bam_ch

    script:
    """
    samtools view -bS ${sam_file} > ${sam_file.baseName}.bam 
    """
}

/*
 * sort bam files by name using samtools
 */
process sort_bam_file {

    cpus = 2
    tag "$state_replicate"
    container "noelnamai/asimov:1.0"

    input:
    set state_replicate, file(bam_file) from aligned_bam_ch

    output:
    set state_replicate, file("${bam_file.baseName}.sorted.bam") into aligned_sorted_bam_ch

    script:
    """
    samtools sort ${bam_file} -n --threads ${task.cpus} -o ${bam_file.baseName}.sorted.bam
    """
}

/*
 * count how many reads map to each gene using htseq-count
 */
process generate_raw_counts {

    cpus = 2
    tag "$state_replicate"
    container "noelnamai/asimov:1.0"

    input:
    file(gtf) from gtf_2_ch
    set state_replicate, file(bam_file) from aligned_sorted_bam_ch

    output:
    file("${state_replicate}.tsv") into raw_counts_ch

    script:
    """
    htseq-count --type exon --idattr gene_id --order name --format bam ${bam_file} ${gtf} > ${state_replicate}.tsv
    """
}

/*
 * detect differentially expressed genes using DESeq2
 */
process detect_dexp_genes {

    cpus = 8
    container "noelnamai/asimov:1.0"

    publishDir "${params.outdir}", mode: "copy", overwrite: true

    input:
    file(count_matrices) from raw_counts_ch.collect()

    output:
    set file("*.png"), file("*.tsv") into dgexp_png_tsv_ch

    script:
    """
    Rscript /opt/dgexp.R --threads ${task.cpus} --results-dir . --alpha 0.1 --log2fc-threshold 1.5
    """
}
