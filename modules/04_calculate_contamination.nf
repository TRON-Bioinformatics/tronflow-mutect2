params.memory_contamination = "16g"
params.cpus_contamination = 2
params.output = 'output'


process CALCULATE_CONTAMINATION {
    cpus params.cpus_contamination
    memory params.memory_contamination
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
    tuple val(name), file(table)

    output:
    tuple val(name), file("${name}.segments.table"), file("${name}.calculatecontamination.table"), emit: contaminationTables

    """
    gatk --java-options '-Xmx${params.memory_contamination}' CalculateContamination \
    --input ${table} \
    -tumor-segmentation ${name}.segments.table \
    --output ${name}.calculatecontamination.table
    """
}