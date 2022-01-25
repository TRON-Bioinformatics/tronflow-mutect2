params.memory_funcotator = "16g"
params.cpus_funcotator = 2
params.db_funcotator = false
params.reference = false
params.reference_version_funcotator = "hg19"
params.output_format_funcotator = "MAF"
params.transcript_selection_mode_funcotator = "CANONICAL"
params.args_funcotator = "--remove-filtered-variants true"

process FUNCOTATOR {

    cpus params.cpus_funcotator
    memory params.memory_funcotator
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)

    input:
    tuple val(name), file(vcf)

    output:
    tuple val(name), val("${params.output}/${name}/${name}.mutect2.funcotated.vcf"), emit: vcf_anno
    file "${name}.mutect2.funcotated.maf"

    """
    gatk --java-options '-Xmx${params.memory_funcotator}' Funcotator \
     --variant ${vcf} \
     --reference ${params.reference} \
     --ref-version ${params.reference_version_funcotator} \
     --data-sources-path ${params.db_funcotator} \
     --output ${name}.mutect2.funcotated.maf \
     --output-file-format ${params.output_format_funcotator} \
     --transcript-selection-mode ${params.transcript_selection_mode_funcotator} \
     ${params.args_funcotator}
    """
}