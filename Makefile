
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
	nextflow main.nf -profile test,conda --output output/test1
	nextflow main.nf -profile test,conda --disable_common_germline_filter --output output/test2
	nextflow main.nf -profile test,conda --input_files test_data/test_input_with_replicates.txt --output output/test3


check:
	test -s output/test1/sample_name/sample_name.mutect2.vcf || { echo "Missing test 1 output file!"; exit 1; }
	test -s output/test2/sample_name/sample_name.mutect2.vcf || { echo "Missing test 2 output file!"; exit 1; }
	test -s output/test3/sample_name_with_replicates/sample_name_with_replicates.mutect2.vcf || { echo "Missing test 3 output file!"; exit 1; }
