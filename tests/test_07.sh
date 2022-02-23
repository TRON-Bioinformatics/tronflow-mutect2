#!/bin/bash

# Test Funcotator functionallity

source bin/assert.sh
output=output/test5

echo -e "sample_name\t"`pwd`"/test_data/SRR8244887.preprocessed.downsampled.bam\t"`pwd`"/test_data/SRR8244836.preprocessed.downsampled.bam" > test_data/test_input.txt
nextflow main.nf -profile test,conda --output $output --input_files test_data/test_input.txt \
  --reference_version_funcotator hg19 \
  --funcotator /home/priesgo/funcotator/funcotator_dataSources.v1.7.20200521s

test -s $output/sample_name/sample_name.mutect2.funcotated.maf || { echo "Missing output VCF file!"; exit 1; }