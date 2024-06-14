#!/bin/bash


source bin/assert.sh
output=output/test10

nextflow main.nf -profile test,conda,ci --output $output \
--input_name sample_name \
--input_tumor_bam `pwd`/test_data/SRR8244887.preprocessed.downsampled.bam \
--input_normal_bam `pwd`/test_data/SRR8244836.preprocessed.downsampled.bam

test -s $output/sample_name/sample_name.mutect2.vcf || { echo "Missing output VCF file!"; exit 1; }