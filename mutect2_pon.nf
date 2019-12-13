#!/usr/bin/env nextflow

gatk4_jar = "/code/gatk/4.1.3.0/gatk-package-4.1.3.0-local.jar"
// this version is needed to build a PON from multiple VCFs, it may change in the future
gatk40_jar = "/code/gatk/4.0.12.0/gatk-package-4.0.12.0-local.jar"
picard_jar = "/code/picard/2.21.2/picard.jar"

params.help= false
params.input_files = false
params.reference = "/projects/data/gatk_bundle/hg19/ucsc.hg19.fasta"							// TODO: remove this hard coded bit
params.intervals = "/projects/data/gatk_bundle/hg19/hg19_refseq_exons.sorted.merged.bed"
params.output = 'output'

def helpMessage() {
    log.info"""
Usage:
    mutect2_pon.nf --input_files input_files

This workflow aims to compute a panel of normals to be used with MuTect2

Input:
    * input_files: the path to a file containing in each row the sample name and the path to a BAM file to be included in the PON
    	example:
	sample1	/path/to/sample1.bam
	sample2	/path/to/sample2.bam
	NOTE: the sample name must be set in the @SN annotation

Optional input:
    * output: the folder where to publish output

Output:
    * Output combined VCF pon.vcf
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'bam'], sep: "\t")
    .map { row-> tuple(row.name, file(row.bam), file(row.bam.replaceAll(/.bam$/, ".bai"))) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}

if (params.intervals) {
  Channel
    .fromPath(params.intervals)
    .splitText(by: 30000, file: true)
    .set { intervals }
}

process mutect2Pon {
    cpus 2
    memory '16g'
    module 'java/1.8.0'
    errorStrategy 'finish'

    input:
    	set val(name), file(bam), file(bai), file(interval) from input_files.combine(intervals)

    output:
	    set val("${bam.baseName}"), file("${bam.baseName}.${interval.baseName}.mutect.vcf")  into mutect_vcfs

    """
    mkdir -p `pwd`/scratch/tmp
    java -Xmx16g -Djava.io.tmpdir=`pwd`/scratch/tmp -jar ${gatk4_jar} \
    Mutect2 \
    --reference ${params.reference} \
    --intervals ${interval} \
    --input ${bam} \
    --tumor-sample ${name} \
    --max-mnp-distance 0 \
    --output ${bam.baseName}.${interval.baseName}.mutect.vcf
    """
}

process gatherVcfs {
    cpus 1
    memory '32g'
    module 'java/1.8.0'
    publishDir "${params.output}", mode: "copy"

    input:
    	set name, file(vcf_list) from mutect_vcfs.groupTuple()	// group by name

    output:
    	file("${name}.vcf") into whole_vcfs

    script:
      	// NOTE: the input VCFs need to be provided in order by genomic coordinates
    	input_vcfs = "$vcf_list".split(" ")
        .sort{ a, b -> a.tokenize(".")[-3].toInteger().compareTo b.tokenize(".")[-3].toInteger() }
        .collect{"INPUT=" + it}.join(" ")
    	"""
      	mkdir -p `pwd`/scratch/tmp
    	java -Xmx32g -Djava.io.tmpdir=`pwd`/scratch/tmp  -jar $picard_jar \
    	GatherVcfs \
      	${input_vcfs} \
      	OUTPUT=${name}.vcf
	    """
}

process createPON {
    cpus 1
    memory '32g'
    module 'java/1.8.0'
    publishDir "${params.output}", mode: "copy"

    input:
    	file(vcf_list) from whole_vcfs.collect()

    output:
    	file("pon.vcf")
      file("pon.vcf.idx")

    script:
    	input_vcfs = "$vcf_list".split(" ").collect{"--vcfs " + it}.join(" ")
    	"""
      # combines VCFs and keeps variants occuring in at least two VCFs
	    mkdir -p `pwd`/scratch/tmp
    	java -Xmx32g -Djava.io.tmpdir=`pwd`/scratch/tmp  -jar $gatk40_jar \
    	CreateSomaticPanelOfNormals ${input_vcfs} \
      --output pon.vcf
    	"""
}
