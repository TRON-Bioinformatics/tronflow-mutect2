
all : clean test check

clean:
	rm -rf output
	#rm -rf work
	rm -f report.html*
	rm -f timeline.html*
	rm -f trace.txt*
	rm -f dag.dot*
	rm -f .nextflow.log*
	rm -rf .nextflow*


test:
	echo "sample_name\t"`pwd`"/test_data/TESTX_S1_L001.bam\t"`pwd`"/test_data/TESTX_S1_L002.bam" > test_data/test_input.txt
	nextflow main.nf -profile test,conda --output output/test1 --input_files test_data/test_input.txt
	nextflow main.nf -profile test,conda --disable_common_germline_filter --output output/test2 --input_files test_data/test_input.txt
	echo "sample_name_with_replicates\t"`pwd`"/test_data/TESTX_S1_L001.bam,"`pwd`"/test_data/TESTX_S1_L001.bam\t"`pwd`"/test_data/TESTX_S1_L002.bam,"`pwd`"/test_data/TESTX_S1_L002.bam" > test_data/test_input_with_replicates.txt
	nextflow main.nf -profile test,conda --input_files test_data/test_input_with_replicates.txt --output output/test3
	nextflow main.nf -profile test,conda --output output/test4 --input_files test_data/test_input.txt --intervals false


check:
	test -s output/test1/sample_name/sample_name.mutect2.vcf || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/sample_name/sample_name.mutect2.vcf || { echo "Missing test 2 output file!"; exit 1; }
	test -s output/test3/sample_name_with_replicates/sample_name_with_replicates.mutect2.vcf || { echo "Missing test 3 output file!"; exit 1; }
	test -s output/test4/sample_name/sample_name.mutect2.vcf || { echo "Missing test 4 output file!"; exit 1; }

