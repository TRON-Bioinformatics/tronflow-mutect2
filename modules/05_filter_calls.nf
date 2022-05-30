params.memory_filter = "16g"
params.output = 'output'
params.reference = false
params.args_filter = ""


process FILTER_CALLS {
    cpus 2
    memory params.memory_filter
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)

    input:
    tuple val(name), path(segments_table), path(contamination_table), path(model), path(unfiltered_vcf), path(vcf_stats)

    output:
    tuple val(name), val("${params.output}/${name}/${name}.mutect2.vcf"), emit: final_vcfs
    tuple val(name), path("${name}.mutect2.vcf"), emit: anno_input
    path "${name}.mutect2.vcf"

    script:
    segments_table_param = segments_table.checkIfExists() ? "--tumor-segmentation ${segments_table}" : ""
    contamination_table_param = contamination_table.checkIfExists() ? "--contamination-table ${contamination_table}" : ""
    """
    gatk --java-options '-Xmx${params.memory_filter}' FilterMutectCalls \
    -V ${unfiltered_vcf} \
    --reference ${params.reference} \
    ${segments_table_param} \
    ${contamination_table_param} \
    --ob-priors ${model} \
    --output ${name}.mutect2.vcf ${params.args_filter}
    """
}