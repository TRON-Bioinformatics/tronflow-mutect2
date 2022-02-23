params.memory_pileup = "32g"
params.output = 'output'
params.gnomad = false


process PILEUP_SUMMARIES {
    cpus 2
    memory params.memory_pileup
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)

    input:
    tuple val(name), val(tumor_bam), val(normal_bam)

    output:
    tuple val("${name}"), file("${name}.pileupsummaries.table"), emit: pileupsummaries

    script:
    tumor_inputs = tumor_bam.split(",").collect({v -> "--input $v"}).join(" ")
    """
    gatk --java-options '-Xmx${params.memory_pileup}' GetPileupSummaries  \
    --intervals ${params.gnomad} \
    --variant ${params.gnomad} \
    ${tumor_inputs} \
    --output ${name}.pileupsummaries.table
    """
}