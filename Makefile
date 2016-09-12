all:
	gcc quantum_sim.c -o quantum_sim -lm -g -std=c99

clean:
	rm -f *.o
	rm -f quantum_sim
