params.memory_mutect2 = "16g"
params.cpus_mutect2 = 2
params.output = 'output'
params.gnomad = false
params.pon = false
params.disable_common_germline_filter = false
params.reference = false
params.intervals = false


process MUTECT2 {
    cpus params.cpus_mutect2
    memory params.memory_mutect2
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)

    input:
    tuple val(name), val(tumor_bam), val(normal_bam)

    output:
    tuple val("${name}"), file("${name}.mutect2.unfiltered.vcf"), file("${name}.mutect2.unfiltered.vcf.stats"), emit: unfiltered_vcfs
    tuple val("${name}"), file("${name}.f1r2.tar.gz"), emit: f1r2_stats

    script:
    normal_panel_option = params.pon ? "--panel-of-normals ${params.pon}" : ""
    germline_filter = params.disable_common_germline_filter ? "" : "--germline-resource ${params.gnomad}"
    normal_inputs = normal_bam.split(",").collect({v -> "--input $v"}).join(" ")
    tumor_inputs = tumor_bam.split(",").collect({v -> "--input $v"}).join(" ")
    intervals_option = params.intervals ? "--intervals ${params.intervals}" : ""
    """
    normalRGSM="\$(samtools view -H ${normal_bam} | grep -oP '(?<=SM:)[^ |\t]*')"
    tumorRGSM="\$(samtools view -H ${tumor_bam} | grep -oP '(?<=SM:)[^ |\t]*')"
  
    gatk --java-options '-Xmx${params.memory_mutect2}' Mutect2 \
    --reference ${params.reference} \
    ${intervals_option} \
    ${germline_filter} \
    ${normal_panel_option} \
    ${normal_inputs} --normal-sample \${normalRGSM} \
    ${tumor_inputs} --tumor-sample \${tumorRGSM} \
    --output ${name}.mutect2.unfiltered.vcf \
    --f1r2-tar-gz ${name}.f1r2.tar.gz
    """
}