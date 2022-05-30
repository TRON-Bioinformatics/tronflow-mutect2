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
params.input_name = false
params.input_tumor_bam = false
params.input_normal_bam = false
params.reference = false
params.gnomad = false
params.output = 'output'
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

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'tumor_bam', 'normal_bam'], sep: "\t")
    .map{ row-> tuple(row.name, row.tumor_bam, row.normal_bam) }
    .set { input_files }
} else if (params.input_name && params.input_tumor_bam && params.input_normal_bam) {
   Channel
    .fromList([tuple(params.input_name, params.input_tumor_bam, params.input_normal_bam)])
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}

workflow {

    MUTECT2(input_files)
    LEARN_READ_ORIENTATION_MODEL(MUTECT2.out.f1r2_stats)

    if (params.gnomad) {
        PILEUP_SUMMARIES(input_files)
        CALCULATE_CONTAMINATION(PILEUP_SUMMARIES.out.pileupsummaries)
        FILTER_CALLS(
            CALCULATE_CONTAMINATION.out.contaminationTables.join(
                LEARN_READ_ORIENTATION_MODEL.out.read_orientation_model).join(MUTECT2.out.unfiltered_vcfs))
    }
    else {
        FILTER_CALLS(
            input_files.map{ row-> tuple(row[0], file("dummy"), file("dummy2")) }.join(
                LEARN_READ_ORIENTATION_MODEL.out.read_orientation_model).join(MUTECT2.out.unfiltered_vcfs))
    }

    FILTER_CALLS.out.final_vcfs.map {it.join("\t")}.collectFile(name: "${params.output}/mutect2_output_files.txt", newLine: true)
    if(params.funcotator){
        FUNCOTATOR(FILTER_CALLS.out.anno_input)
    }
}
