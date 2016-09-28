all:
	gcc quantum_sim.c -o quantum_sim -lm -std=c99

openmp:
	gcc quantum_sim.c -o quantum_sim_openmp -lm -std=c99 -fopenmp

phi:
	icc -qopenmp -O3 -qopt-report=3 -qopt-report-phase=vec -mmic quantum_sim.c -o quantum_sim_phi -lm -std=c99 -ldl -lpthread

phi_no_openmp:
	icc -O3 -qopt-report=3 -qopt-report-phase=vec -mmic quantum_sim.c -o quantum_sim_phi_no_openmp -lm -std=c99 -ldl -lpthread

clean:
	rm -f *.o
	rm -f quantum_sim quantum_sim_openmp quantum_sim_phi quantum_sim_phi_no_openmp
