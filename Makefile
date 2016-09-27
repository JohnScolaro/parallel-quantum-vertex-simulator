all:
	gcc quantum_sim.c -o quantum_sim -lm -std=c99

openmp:
	gcc quantum_sim.c -o quantum_sim -lm -std=c99 -fopenmp

clean:
	rm -f *.o
	rm -f quantum_sim
