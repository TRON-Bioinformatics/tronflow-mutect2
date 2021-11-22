#!/bin/bash


source tests/assert.sh
output=output/test2

nextflow main.nf -profile test,conda --disable_common_germline_filter --output $output --input_files test_data/test_input.txt

test -s $output/sample_name/sample_name.mutect2.vcf || { echo "Missing output VCF file!"; exit 1; }