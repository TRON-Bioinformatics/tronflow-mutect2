#!/usr/bin/env nextflow


params.help= false
params.input_files = false
params.reference = false
params.intervals = false
params.gnomad = false
params.output = 'output'
params.pon = false
params.memory_mutect2 = "16g"
params.cpus_mutect2 = 2
params.memory_read_orientation = "16g"
params.cpus_read_orientation = 2
params.memory_pileup = "32g"
params.cpus_pileup = 2
params.memory_contamination = "16g"
params.cpus_contamination = 2
params.memory_filter = "16g"
params.cpus_filter = 2
params.disable_common_germline_filter = false

def helpMessage() {
    log.info params.help_message
}

if (params.help) {
    helpMessage()
    exit 0
}
if (!params.reference) {
    log.error "--reference is required"
    exit 1
}
if (!params.gnomad) {
    log.error "--gnomad is required"
    exit 1
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'tumor_bam', 'normal_bam'], sep: "\t")
    .map{ row-> tuple(row.name, row.tumor_bam, row.normal_bam) }
    .set { input_files }

  Channel
      .fromPath(params.input_files)
      .splitCsv(header: ['name', 'tumor_bam', 'normal_bam'], sep: "\t")
      .map{ row-> tuple(row.name, row.tumor_bam) }
      .set { tumor_bams }
} else {
  exit 1, "Input file not specified!"
}

process mutect2 {
    cpus params.cpus_mutect2
    memory params.memory_mutect2
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
    	set name, tumor_bam, normal_bam from input_files

    output:
	    set val("${name}"), file("${name}.mutect2.unfiltered.vcf"), file("${name}.mutect2.unfiltered.vcf.stats") into unfiltered_vcfs
        set val("${name}"), file("${name}.f1r2.tar.gz") into f1r2_stats

    script:
    	normal_panel_option = params.pon ? "--panel-of-normals ${params.pon}" : ""
    	germline_filter = params.disable_common_germline_filter ? "" : "--germline-resource ${params.gnomad}"
    	normal_inputs = normal_bam.split(",").collect({v -> "--input $v"}).join(" ")
    	tumor_inputs = tumor_bam.split(",").collect({v -> "--input $v"}).join(" ")
    	intervals_option = params.intervals ? "--intervals ${params.intervals}" : ""
	"""
    gatk --java-options '-Xmx${params.memory_mutect2}' Mutect2 \
	--reference ${params.reference} \
	${intervals_option} \
	${germline_filter} \
	${normal_panel_option} \
	${normal_inputs} --normal-sample normal \
    ${tumor_inputs} --tumor-sample tumor \
	--output ${name}.mutect2.unfiltered.vcf \
    --f1r2-tar-gz ${name}.f1r2.tar.gz
	"""
}

process learnReadOrientationModel {
  cpus params.cpus_read_orientation
  memory params.memory_read_orientation
  tag "${name}"
  publishDir "${params.output}/${name}", mode: "copy"

  input:
    set name, file(f1r2_stats) from f1r2_stats

  output:
    set name, file("${name}.read-orientation-model.tar.gz") into read_orientation_model

  """
  gatk --java-options '-Xmx${params.memory_read_orientation}' LearnReadOrientationModel \
  --input ${f1r2_stats} \
  --output ${name}.read-orientation-model.tar.gz
  """
}

process pileUpSummaries {
    cpus params.cpus_pileup
    memory params.memory_pileup
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
    	set name, tumor_bam from tumor_bams

    output:
	   set val("${name}"), file("${name}.pileupsummaries.table") into pileupsummaries

    script:
    tumor_inputs = tumor_bam.split(",").collect({v -> "--input $v"}).join(" ")
    intervals_option = params.intervals ? "--intervals ${params.intervals}" : ""
	"""
    gatk --java-options '-Xmx${params.memory_pileup}' GetPileupSummaries  \
	${intervals_option} \
	--variant ${params.gnomad} \
	${tumor_inputs} \
	--output ${name}.pileupsummaries.table
	"""
}

process calculateContamination {
    cpus params.cpus_contamination
    memory params.memory_contamination
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
      set name, file(table) from pileupsummaries

    output:
      set name, file("${name}.segments.table"), file("${name}.calculatecontamination.table") into contaminationTables

    """
    gatk --java-options '-Xmx${params.memory_contamination}' CalculateContamination \
    --input ${table} \
    -tumor-segmentation ${name}.segments.table \
    --output ${name}.calculatecontamination.table
    """
}

process filterCalls {
    cpus params.cpus_filter
    memory params.memory_filter
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
      set name, file(segments_table), file(contamination_table), file(model),
      file(unfiltered_vcf), file(vcf_stats) from contaminationTables.join(read_orientation_model).join(unfiltered_vcfs)

    output:
      set name, val("${params.output}/${name}/${name}.mutect2.vcf") into final_vcfs
      file "${name}.mutect2.vcf"

    """
    gatk --java-options '-Xmx${params.memory_filter}' FilterMutectCalls \
    -V ${unfiltered_vcf} \
    --reference ${params.reference} \
    --tumor-segmentation ${segments_table} \
    --contamination-table ${contamination_table} \
    --ob-priors ${model} \
    --output ${name}.mutect2.vcf
    """
}

final_vcfs.map {it.join("\t")}.collectFile(name: "${params.output}/mutect2_output_files.txt", newLine: true)
