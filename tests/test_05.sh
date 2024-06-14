#!/bin/bash
##################################################################
# Error condition: RGSM tag is the same in tumor and normal   ####
##################################################################

source bin/assert.sh
output=output/test5

echo -e "sample_name\t"`pwd`"/test_data/SRR8244887.preprocessed.downsampled.bam\t"`pwd`"/test_data/SRR8244836.preprocessed.downsampled.bam" > test_data/test_input.txt

{ # try
    nextflow main.nf -profile test,conda,ci --output $output --input_files test_data/test_input.txt &&
    assert_true false "Error condition not captured"
} || { # catch
    assert_true true
}