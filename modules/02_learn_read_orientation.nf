params.memory_read_orientation = "16g"
params.output = 'output'


process LEARN_READ_ORIENTATION_MODEL {
  cpus 2
  memory params.memory_read_orientation
  tag "${name}"
  publishDir "${params.output}/${name}", mode: "copy"

  conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)

  input:
  tuple val(name), file(f1r2_stats)

  output:
  tuple val(name), file("${name}.read-orientation-model.tar.gz"), emit: read_orientation_model

  """
  gatk --java-options '-Xmx${params.memory_read_orientation}' LearnReadOrientationModel \
  --input ${f1r2_stats} \
  --output ${name}.read-orientation-model.tar.gz
  """
}
