#!/usr/bin/env nextflow


params.help= false
params.input_files = false
params.reference = "/projects/data/gatk_bundle/hg19/ucsc.hg19.fasta"
// this BED does not have those contigs that were removed from the reference when aligning
params.intervals = "/projects/data/gatk_bundle/hg19/hg19_refseq_exons.sorted.for_strelka.bed.gz"
params.gnomad = "/projects/data/gatk_bundle/hg19/gnomad.exomes.r2.1.1.sites.PASS.only_af.vcf.bgz"
params.output = 'output'
params.pon = false
params.memory = "16g"
params.cpus = 2

def helpMessage() {
    log.info"""
Usage:
    nextflow run tron-bioinformatics/tronflow-mutect2 -profile conda --input_files input_files [--reference reference.fasta]

This workflow is based on the implementation at /code/iCaM/scripts/mutect2_ID.sh

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name, tumor bam and normal bam
    The input file does not have header!
    Example input file:
    name1	tumor_bam1	normal_bam1
    name2	tumor_bam2	normal_bam2

Optional input:
    * reference: path to the FASTA genome reference (indexes expected *.fai, *.dict)
    * intervals: path to a BED file containing the regions to analyse
    * gnomad: path to the gnomad VCF
    * NOTE: if any of the above parameters is not provided, default hg19 resources will be used
    * output: the folder where to publish output
    * memory: the ammount of memory used by each job (default: 16g)
    * cpus: the number of CPUs used by each job (default: 2)

Output:
    * Output VCF
    * Other intermediate files
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
    .splitCsv(header: ['name', 'tumor_bam', 'normal_bam'], sep: "\t")
    .map{ row->
	    tuple(row.name,
	    file(row.tumor_bam), file(row.tumor_bam + ".bai"),
	    file(row.normal_bam), file(row.normal_bam + ".bai")
	    ) }
    .set { input_files }

  Channel
      .fromPath(params.input_files)
      .splitCsv(header: ['name', 'tumor_bam', 'normal_bam'], sep: "\t")
      .map{ row->
  	    tuple(row.name,
  	    file(row.tumor_bam), file(row.tumor_bam + ".bai")) }
      .set { tumor_bams }
} else {
  exit 1, "Input file not specified!"
}

process mutect2 {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
    	set name, file(tumor_bam), file(tumor_bai), file(normal_bam), file(normal_bai) from input_files

    output:
	    set val("${name}"), file("${name}.unfiltered.vcf"), file("${name}.unfiltered.vcf.stats") into unfiltered_vcfs
      set val("${name}"), file("${name}.f1r2.tar.gz") into f1r2_stats

    script:
    	normal_panel_option = params.pon ? "--panel-of-normals ${params.pon}" : ""
	"""
    gatk --java-options '-Xmx${params.memory}' Mutect2 \
	--reference ${params.reference} \
	--intervals ${params.intervals} \
	--germline-resource ${params.gnomad} \
	${normal_panel_option} \
	--input ${normal_bam} \
  --normal-sample normal \
  --input ${tumor_bam} \
  --tumor-sample tumor \
	--output ${name}.unfiltered.vcf \
  --f1r2-tar-gz ${name}.f1r2.tar.gz
	"""
}

process learnReadOrientationModel {
  cpus params.cpus
  memory params.memory
  tag "${name}"
  publishDir "${params.output}", mode: "copy"

  input:
    set name, file(f1r2_stats) from f1r2_stats

  output:
    set name, file("${name}.read-orientation-model.tar.gz") into read_orientation_model

  """
  gatk --java-options '-Xmx${params.memory}' LearnReadOrientationModel \
  --input ${f1r2_stats} \
  --output ${name}.read-orientation-model.tar.gz
  """
}

process pileUpSummaries {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
    	set name, file(tumor_bam), file(tumor_bai) from tumor_bams

    output:
	   set val("${name}"), file("${name}.pileupsummaries.table") into pileupsummaries

    script:
	"""
    gatk --java-options '-Xmx${params.memory}' GetPileupSummaries  \
	--intervals ${params.gnomad} \
	--variant ${params.gnomad} \
	--input ${tumor_bam} \
	--output ${name}.pileupsummaries.table
	"""
}

process calculateContamination {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
      set name, file(table) from pileupsummaries

    output:
      set name, file("${name}.segments.table"), file("${name}.calculatecontamination.table") into contaminationTables

    """
    gatk --java-options '-Xmx${params.memory}' CalculateContamination \
    --input ${table} \
    -tumor-segmentation ${name}.segments.table \
    --output ${name}.calculatecontamination.table
    """
}

process filterCalls {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
      set name, file(segments_table), file(contamination_table), file(model),
      file(unfiltered_vcf), file(vcf_stats) from contaminationTables.join(read_orientation_model).join(unfiltered_vcfs)

    output:
      set name, file("${name}.vcf") into final_vcfs

    """
    gatk --java-options '-Xmx${params.memory}' FilterMutectCalls \
    -V ${unfiltered_vcf} \
    --reference ${params.reference} \
    --tumor-segmentation ${segments_table} \
    --contamination-table ${contamination_table} \
    --ob-priors ${model} \
    --output ${name}.vcf
    """
}

final_vcfs.map {it.join("\t")}.collectFile(name: "${params.output}/mutect2_output_files.txt", newLine: true)
