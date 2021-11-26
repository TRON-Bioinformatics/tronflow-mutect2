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

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0 bioconda::samtools=1.12" : null)

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
    normalRGSMs = normal_bam.split(",").collect({v -> "\$(samtools view -H $v | grep -oP '(?<=SM:)[^ |\\t]*')"})
    normalRGSM = normalRGSMs.first()
    tumorRGSMs = tumor_bam.split(",").collect({v -> "\$(samtools view -H $v | grep -oP '(?<=SM:)[^ |\\t]*')"})
    tumorRGSM = tumorRGSMs.first()
    intervals_option = params.intervals ? "--intervals ${params.intervals}" : ""
    """
    # sanity checks on the RGSM
    source assert.sh

    assert_eq \$(echo ${tumorRGSMs} | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/, /\\n/g' | sort | uniq | wc -l) 1 "All tumor BAMs RGSM tags must be equal"
    assert_eq \$(echo ${normalRGSMs} | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/, /\\n/g' | sort | uniq | wc -l) 1 "All normal BAMs RGSM tags must be equal"
    assert_not_eq "${normalRGSM}" "${tumorRGSM}" "Tumor and normal RGSM must be different!"

    gatk --java-options '-Xmx${params.memory_mutect2}' Mutect2 \
    --reference ${params.reference} \
    ${intervals_option} \
    ${germline_filter} \
    ${normal_panel_option} \
    ${normal_inputs} --normal-sample ${normalRGSM} \
    ${tumor_inputs} --tumor-sample ${tumorRGSM} \
    --output ${name}.mutect2.unfiltered.vcf \
    --f1r2-tar-gz ${name}.f1r2.tar.gz
    """
}
