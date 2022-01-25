#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MUTECT2 } from './modules/01_mutect2'
include { LEARN_READ_ORIENTATION_MODEL } from './modules/02_learn_read_orientation'
include { PILEUP_SUMMARIES } from './modules/03_pileup_summary'
include { CALCULATE_CONTAMINATION } from './modules/04_calculate_contamination'
include { FILTER_CALLS } from './modules/05_filter_calls'
include { FUNCOTATOR } from './modules/06_annotate'

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
params.funcotator = false

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
if (!params.db_funcotator & params.funcotator) {
    log.error "--db_funcotator is required"
    exit 1
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'tumor_bam', 'normal_bam'], sep: "\t")
    .map{ row-> tuple(row.name, row.tumor_bam, row.normal_bam) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}

workflow {
    MUTECT2(input_files)
    PILEUP_SUMMARIES(input_files)
    LEARN_READ_ORIENTATION_MODEL(MUTECT2.out.f1r2_stats)
    CALCULATE_CONTAMINATION(PILEUP_SUMMARIES.out.pileupsummaries)
    FILTER_CALLS(
        CALCULATE_CONTAMINATION.out.contaminationTables.join(
            LEARN_READ_ORIENTATION_MODEL.out.read_orientation_model).join(MUTECT2.out.unfiltered_vcfs))

    FILTER_CALLS.out.final_vcfs.map {it.join("\t")}.collectFile(name: "${params.output}/mutect2_output_files.txt", newLine: true)
    if(params.funcotator){ // VCF output param
        FUNCOTATOR(FILTER_CALLS.out.anno_input)
    }

}
