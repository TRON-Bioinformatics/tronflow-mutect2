params.memory_filter = "16g"
params.cpus_filter = 2
params.output = 'output'
params.reference = false
params.args_filter = "--min-allele-fraction 0.05 --min-reads-per-strand 1 --unique-alt-read-count 4" // add to documentation


process FILTER_CALLS {
    cpus params.cpus_filter
    memory params.memory_filter
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)

    input:
    tuple val(name), file(segments_table), file(contamination_table), file(model), file(unfiltered_vcf), file(vcf_stats)

    output:
    tuple val(name), val("${params.output}/${name}/${name}.mutect2.vcf"), emit: final_vcfs
    tuple val(name), file("${name}.mutect2.vcf"), emit: anno_input
    file "${name}.mutect2.vcf"

    """
    gatk --java-options '-Xmx${params.memory_filter}' FilterMutectCalls \
    -V ${unfiltered_vcf} \
    --reference ${params.reference} \
    --tumor-segmentation ${segments_table} \
    --contamination-table ${contamination_table} \
    --ob-priors ${model} \
    --output ${name}.mutect2.vcf ${params.args_filter}
    """
}