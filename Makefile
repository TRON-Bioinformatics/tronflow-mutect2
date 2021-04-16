clean:
	rm -rf output
	rm -rf work
	rm -f report.html*
	rm -f timeline.html*
	rm -f trace.txt*
	rm -f dag.dot*
	rm -f .nextflow.log*
	rm -rf .nextflow*


test:
	nextflow main.nf -profile test,conda
	nextflow main.nf -profile test,conda --disable_common_germline_filter
