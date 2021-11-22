params.memory_read_orientation = "16g"
params.cpus_read_orientation = 2
params.output = 'output'


process LEARN_READ_ORIENTATION_MODEL {
  cpus params.cpus_read_orientation
  memory params.memory_read_orientation
  tag "${name}"
  publishDir "${params.output}/${name}", mode: "copy"

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
