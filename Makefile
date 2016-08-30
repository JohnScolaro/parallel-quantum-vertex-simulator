all: 
	icc -fopenmp -mmic -vec-report=3 -O3 helloflops2.c -o helloflops2
