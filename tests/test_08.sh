#!/bin/bash


source bin/assert.sh
output=output/test8

echo -e "sample_name\t"`pwd`"/test_data/SRR8244887.preprocessed.downsampled.bam\t"`pwd`"/test_data/SRR8244836.preprocessed.downsampled.bam" > test_data/test_input.txt
nextflow main.nf -profile test,conda --output $output --input_files test_data/test_input.txt --enable_bam_output

test -s $output/sample_name/sample_name.mutect2.vcf || { echo "Missing output VCF file!"; exit 1; }
test -s $output/sample_name/sample_name.mutect2.assembled_haplotypes.bam || { echo "Missing output BAM file!"; exit 1; }
test -s $output/sample_name/sample_name.mutect2.assembled_haplotypes.bai || { echo "Missing output BAI file!"; exit 1; }