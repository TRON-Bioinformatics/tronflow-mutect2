params.memory_contamination = "16g"
params.output = 'output'


process CALCULATE_CONTAMINATION {
    cpus 2
    memory params.memory_contamination
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)

    input:
    tuple val(name), path(table)

    output:
    tuple val(name), path("${name}.segments.table"), path("${name}.calculatecontamination.table"), emit: contaminationTables

    """
    gatk --java-options '-Xmx${params.memory_contamination}' CalculateContamination \
    --input ${table} \
    -tumor-segmentation ${name}.segments.table \
    --output ${name}.calculatecontamination.table
    """
}