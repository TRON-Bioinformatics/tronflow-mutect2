#!/bin/bash


source bin/assert.sh
output=output/test9

echo -e "sample_name\t"`pwd`"/test_data/SRR8244887.preprocessed.downsampled.bam\t"`pwd`"/test_data/SRR8244836.preprocessed.downsampled.bam" > test_data/test_input.txt
nextflow main.nf -profile test,conda,ci --output $output --input_files test_data/test_input.txt --gnomad false  --args_filter "--contamination-estimate 0.2"

test -s $output/sample_name/sample_name.mutect2.vcf || { echo "Missing output VCF file!"; exit 1; }