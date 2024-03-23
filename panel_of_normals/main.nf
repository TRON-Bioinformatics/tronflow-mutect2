#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// params definitions section

params.help= false
params.gnomad = false
params.input_files = false
params.intervals = false
params.max_mnp_distance = 0
params.memory_mutect2 = "16g"
params.memory_create_genomics_db = "16g"
params.memory_create_pon = "32g"
params.reference = false
params.output = "output"

// checks of required input arguments

def helpMessage() {
    log.info params.help_message
}

if (params.help) {
    helpMessage()
    exit 0
}

if (params.input_files) {
    Channel
        .fromPath(params.input_files)
        .splitCsv(header: ["name", "bam"], sep: "\t")
        .map { row-> tuple(row.name, row.bam) }
        .set { input_files }
} else {
    exit 1, "The following argument is required: --input_files!"
}

if (! params.reference) {
    exit 1, "The following argument is required: --reference!"
}

if (! params.gnomad) {
    exit 1, "The following argument is required: --gnomad!"
}

// process definitions section

process MUTECT2 {
    cpus 2
    memory params.memory_mutect2
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda "bioconda::gatk4=4.2.6.1 bioconda::samtools=1.12"

    input:
    	tuple val(name), val(bam)

    output:
	    tuple val("${name}"), path("${name}.mutect2.unfiltered.vcf"), emit: vcf
        tuple val("${name}"), path("${name}.mutect2.unfiltered.vcf.idx"), emit: vcf_idx
        tuple val("${name}"), path("${name}.mutect2.unfiltered.vcf.stats"), emit: vcf_stats

    script:
        intervals_option = params.intervals ? "--intervals ${params.intervals}" : ""

        """
        mkdir -p `pwd`/scratch/tmp
        gatk --java-options '-Xmx${params.memory_mutect2}' Mutect2 \\
            --reference ${params.reference} \\
            ${intervals_option} \\
            --input ${bam} \\
            --tumor-sample ${name} \\
            --max-mnp-distance ${params.max_mnp_distance} \\
            --output ${name}.mutect2.unfiltered.vcf \\
            --tmp-dir `pwd`/scratch/tmp
        """
}

process CREATE_GENOMICS_DB {
    cpus 2
    memory params.memory_create_genomics_db
    tag "genomics_db"

    conda "bioconda::gatk4=4.2.6.1 bioconda::samtools=1.12"

    input:
    	val(vcf_list)

    output:
	    path("genomics_db", type: "dir")

    script:
        intervals_option = params.intervals ? "--intervals ${params.intervals}" : ""
        input_vcfs = vcf_list.collect{"--variant " + it}.join(" ")

        """
        mkdir -p `pwd`/scratch/tmp
        gatk --java-options '-Xmx${params.memory_create_genomics_db}' GenomicsDBImport \\
            --reference ${params.reference} \\
            ${intervals_option} \\
            --genomicsdb-workspace-path genomics_db \\
            --tmp-dir `pwd`/scratch/tmp \\
            --genomicsdb-shared-posixfs-optimizations true \\
            --batch-size 50 \\
            --bypass-feature-reader \\
            ${input_vcfs}
        """
}

process CREATE_PON {
    cpus 1
    memory params.memory_create_pon
    tag "pon.vcf"
    publishDir "${params.output}", mode: "copy"

    conda "bioconda::gatk4=4.2.6.1 bioconda::samtools=1.12"

    input:
    	path(genomics_db)

    output:
    	file("pon.vcf")
        file("pon.vcf.idx")

    script:
    	"""
        mkdir -p `pwd`/scratch/tmp
        gatk --java-options '-Xmx${params.memory_create_pon}' CreateSomaticPanelOfNormals \\
            --reference ${params.reference} \\
            --germline-resource ${params.gnomad} \\
            --variant gendb://${genomics_db} \\
            --output pon.vcf \\
            --tmp-dir `pwd`/scratch/tmp
    	"""
}

// workflow definitions section

workflow {
    MUTECT2(input_files)

    // CREATE_GENOMICS_DB(
    //     MUTECT2.out.vcf
    //         .map { name, vcf -> vcf }
    //         .collect()
    // )

    // CREATE_PON(CREATE_GENOMICS_DB.out)
}
