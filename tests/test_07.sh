#!/bin/bash

# Test Funcotator functionallity

source bin/assert.sh
output=output/test5

echo -e "sample_name\t"`pwd`"/test_data/TESTX_S1_L001.bam\t"`pwd`"/test_data/TESTX_S1_L002.bam" > test_data/test_input.txt
nextflow main.nf -profile test,conda --output $output --input_files test_data/test_input.txt --reference_version_funcotator="hg19" --db_funcotator /flash/home/chritzel/sarcoma_analysis/exome/references/funcotator_dataSources.v1.6.20190124s/ --funcotator true

test -s $output/sample_name/sample_name.mutect2.funcotated.maf || { echo "Missing output VCF file!"; exit 1; }