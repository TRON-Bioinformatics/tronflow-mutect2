#!/bin/bash

##############################################################################
# Error condition: RGSM tag is the different betweeen two tumor           ####
# samples or two normal samples                                           ####
##############################################################################

source bin/assert.sh
output=output/test6

echo -e "sample_name\t"`pwd`"/test_data/TESTX_S1_L001.bam,"`pwd`"/test_data/TESTX_S1_L002.bam\t"`pwd`"/test_data/TESTX_S1_L002.bam" > test_data/test_input.txt
{ # try
    nextflow main.nf -profile test,conda --output $output --input_files test_data/test_input.txt &&
    assert_true false "Error condition not captured"
} || { # catch
    assert_true true
}

echo -e "sample_name\t"`pwd`"/test_data/TESTX_S1_L001.bam\t"`pwd`"/test_data/TESTX_S1_L001.bam,"`pwd`"/test_data/TESTX_S1_L002.bam" > test_data/test_input.txt
{ # try
    nextflow main.nf -profile test,conda --output $output --input_files test_data/test_input.txt &&
    assert_true false "Error condition not captured"
} || { # catch
    assert_true true
}