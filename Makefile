
all : clean test

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
	bash tests/test_00.sh
	bash tests/test_01.sh
	bash tests/test_02.sh
	bash tests/test_03.sh
	bash tests/test_04.sh
	bash tests/test_05.sh
	bash tests/test_06.sh
	bash tests/test_08.sh
	bash tests/test_09.sh
	bash tests/test_10.sh
