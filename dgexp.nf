#!/usr/bin/env nextflow

params.help = null

println """\

R N A - S E Q  D I F F E R E N T I A L  E X P R E S I O N 
=========================================================

HISAT2     : 2.1.0
Samtools   : 1.7
Started at : $workflow.start
Results dir: ${params.outdir}
"""
.stripIndent()

/* 
 * parse the input parameters.
 */
genome     = file(params.genome)
annotation = file(params.annot)

/*
 * read fastq.gz files into a channel
 */
Channel.fromPath(params.reads)
    .ifEmpty{error "Cannot find any reads matching: ${params.reads}"}
    .map{file -> [file.name.substring(0, 17), file]}
    .groupTuple(by: 0, sort: true)
    .set{read_pairs_ch}

/*
 * extract exons and splice junctions from the gff3 file
 */
process extract_exons_and_ss {

    container "noelnamai/asimov-dgexp:v1.0"

    cpus = 2

    input:
    file annotation
      
    output:
    set file("${annotation.baseName}.exons.tsv"), file("${annotation.baseName}.splicesites.tsv") into extracted_exons_and_ss_ch

    script:
    """
    gt gff3_to_gtf -o ${annotation.baseName}.gtf ${annotation}
    
    hisat2_extract_exons.py ${annotation.baseName}.gtf > ${annotation.baseName}.exons.tsv
    hisat2_extract_splice_sites.py ${annotation.baseName}.gtf > ${annotation.baseName}.splicesites.tsv
    """
}

/*
 * build the genome index for mapping using hisat2
 */
process build_genome_index {
    
    container "noelnamai/asimov-dgexp:v1.0"

    cpus = 2

    input:
    file genome
    set file(extracted_exons), file(extracted_ss) from extracted_exons_and_ss_ch
      
    output:
    file("${genome}.index.*") into genome_index_ch

    script:
    """
    hisat2-build --exon ${extracted_exons} --ss ${extracted_ss} -p ${task.cpus} ${genome} ${genome}.index
    """
}

/*
 * align reads to the reference genome using hisat2
 */
process map_reads_to_reference {
    
    tag "$state_rep"

    container "noelnamai/asimov-dgexp:v1.0"

    cpus = 2

    input:
    file genome
    file index from genome_index_ch
    set state_rep, file(reads) from read_pairs_ch

    output:
    set state_rep, file("${state_rep}.sam") into aligned_sam_ch

    script:
    """
    hisat2 --dta -p ${task.cpus} -x ${genome}.index -U ${reads[0]},${reads[1]} -S ${state_rep}.sam
    """
}

/*
 * convert sam to sorted bam files using samtools
 */
process convert_sam_to_bam {

    tag "$state_rep"

    container "noelnamai/asimov-dgexp:v1.0"

    cpus = 2

    input:
    set state_rep, file(sam_file) from aligned_sam_ch

    output:
    set state_rep, file("${sam_file.baseName}.bam") into aligned_bam_ch

    script:
    """
    samtools view -bS ${sam_file} > ${sam_file.baseName}.bam 
    """
}

/*
 * sort bam files using samtools
 */
process sort_bam_file {

    tag "$state_rep"

    container "noelnamai/asimov-dgexp:v1.0"

    cpus = 2

    input:
    set state_rep, file(bam_file) from aligned_bam_ch

    output:
    set state_rep, file("${bam_file.baseName}.sorted.bam") into aligned_sorted_bam_ch

    script:
    """
    samtools sort ${bam_file} -n --threads ${task.cpus} -o ${bam_file.baseName}.sorted.bam
    """
}

process generate_raw_counts {

    tag "$state_rep"

    container "noelnamai/asimov-dgexp:v1.0"

    cpus = 2

    input:
    file annotation
    set state_rep, file(bam_file) from aligned_sorted_bam_ch

    output:
    file("${bam_file.baseName}.tsv") into raw_counts_ch

    script:
    """
    gt gff3_to_gtf -o ${annotation.baseName}.gtf ${annotation}

    htseq-count --format bam ${bam_file} ${annotation.baseName}.gtf > ${bam_file.baseName}.tsv
    """
}

raw_counts_ch
    .collect()
    .println()