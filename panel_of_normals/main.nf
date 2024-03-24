#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// params definitions section

params.help= false
params.gnomad = false
params.input_files = false
params.intervals = false
params.max_mnp_distance = 0
params.memory_mutect2 = "16g"
params.memory_gather_vcfs = "32g"
params.memory_genomicsdb_import = "16g"
params.memory_create_pon = "32g"
params.reference = false
params.output = "output"
params.genomicsdb_output = "output"

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

if (params.intervals) {
    Channel
        .fromPath(params.intervals)
        .splitText(by: 30000, file: true)
        .set { intervals }
} else {
    exit 1, "The following argument is required: --intervals!"
}

// process definitions section

process MUTECT2 {
    cpus 2
    memory params.memory_mutect2
    tag "${name}-${interval}"

    conda "bioconda::gatk4=4.2.6.1 bioconda::samtools=1.12"

    input:
    	tuple val(name), val(bam), path(interval)

    output:
	    tuple val("${name}"), path("${name}.${interval.baseName}.mutect2.unfiltered.vcf.gz")

    script:
        """
        mkdir -p `pwd`/scratch/tmp
        gatk --java-options '-Xmx${params.memory_mutect2}' Mutect2 \\
            --reference ${params.reference} \\
            --intervals ${interval} \\
            --input ${bam} \\
            --tumor-sample ${name} \\
            --max-mnp-distance ${params.max_mnp_distance} \\
            --output ${name}.${interval.baseName}.mutect2.unfiltered.vcf.gz \\
            --tmp-dir `pwd`/scratch/tmp
        """
}

process GATHER_VCFS {
    cpus 1
    memory params.memory_gather_vcfs
    tag "${name}"

    conda "bioconda::gatk4=4.2.6.1 bioconda::samtools=1.12"

    input:
    	tuple val(name), val(vcf_list)

    output:
	    tuple val("${name}"), path("${name}.mutect2.unfiltered.vcf.gz"), path("${name}.mutect2.unfiltered.vcf.gz.tbi")

    script:
        // NOTE: the input VCFs need to be provided in order by genomic coordinates
    	input_vcfs = vcf_list
            .sort { a, b -> file(a).baseName.compareTo(file(b).baseName) }
            .collect{"-INPUT " + it}.join(" ")
        """
        mkdir -p `pwd`/scratch/tmp
        gatk --java-options '-Xmx${params.memory_gather_vcfs} -Djava.io.tmpdir=`pwd`/scratch/tmp' GatherVcfs \\
            ${input_vcfs} \\
            -OUTPUT ${name}.mutect2.unfiltered.vcf.gz
        tabix -p vcf ${name}.mutect2.unfiltered.vcf.gz
        """
}

process GENOMICSDB_IMPORT {
    cpus 2
    memory params.memory_genomicsdb_import
    tag "${name}"

    conda "bioconda::gatk4=4.2.6.1 bioconda::samtools=1.12"

    input:
    	tuple val(name), path(vcf), path(vcf_tbi)

    output:
	    val("Import completed")

    script:
        if (! file("${params.genomicsdb_output}/pon.gdb").exists() ||  
            file("${params.genomicsdb_output}/pon.gdb").list().length == 0) {
        """
        mkdir -p `pwd`/scratch/tmp
        gatk --java-options '-Xmx${params.memory_genomicsdb_import}' GenomicsDBImport \\
            --reference ${params.reference} \\
            --intervals ${params.intervals} \\
            --genomicsdb-workspace-path ${params.genomicsdb_output}/pon.gdb \\
            --tmp-dir `pwd`/scratch/tmp \\
            --genomicsdb-shared-posixfs-optimizations true \\
            --batch-size 50 \\
            --bypass-feature-reader \\
            --merge-input-intervals \\
            --variant ${vcf}
        """
        } else {
        """
        mkdir -p `pwd`/scratch/tmp
        gatk --java-options '-Xmx${params.memory_genomicsdb_import}' GenomicsDBImport \\
            --reference ${params.reference} \\
            --intervals ${params.intervals} \\
            --genomicsdb-update-workspace-path ${params.genomicsdb_output}/pon.gdb \\
            --tmp-dir `pwd`/scratch/tmp \\
            --genomicsdb-shared-posixfs-optimizations true \\
            --batch-size 50 \\
            --bypass-feature-reader \\
            --merge-input-intervals \\
            --variant ${vcf}
        """
        }
}

process CREATE_PON {
    cpus 1
    memory params.memory_create_pon
    tag "pon.vcf"
    publishDir "${params.output}", mode: "copy"

    conda "bioconda::gatk4=4.2.6.1 bioconda::samtools=1.12"

    input:
        val(genomicsdb_import_completed)

    output:
    	file("pon.vcf")
        file("pon.vcf.idx")

    script:
    	"""
        mkdir -p `pwd`/scratch/tmp
        gatk --java-options '-Xmx${params.memory_create_pon}' CreateSomaticPanelOfNormals \\
            --reference ${params.reference} \\
            --germline-resource ${params.gnomad} \\
            --variant gendb://${params.genomicsdb_output}/pon.gdb \\
            --output pon.vcf \\
            --tmp-dir `pwd`/scratch/tmp
    	"""
}

// workflow definitions section

workflow {
    MUTECT2(input_files.combine(intervals))

    GATHER_VCFS(MUTECT2.out.groupTuple())

    GENOMICSDB_IMPORT(GATHER_VCFS.out)

    CREATE_PON(GENOMICSDB_IMPORT.out.collect())
}
