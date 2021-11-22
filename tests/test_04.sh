#!/bin/bash


source tests/assert.sh
output=output/test4

nextflow main.nf -profile test,conda --output $output --input_files test_data/test_input.txt --intervals false

test -s $output/sample_name/sample_name.mutect2.vcf || { echo "Missing output VCF file!"; exit 1; }