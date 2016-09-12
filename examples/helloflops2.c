/*
 * An example to try to get lots of flops on the Intel Xeon Phi.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>

double dtime(void) {
	double tseconds = 0.0;
	struct timeval mytime;

	gettimeofday(&mytime, (struct timezone*)0);
	tseconds = (double) (mytime.tv_sec + mytime.tv_usec * 1.0e-6);

	return tseconds;
}

#define FLOPS_ARRAY_SIZE (1024 * 1024)
#define MAXFLOPS_ITERS 100000000
#define LOOP_COUNT 128
#define FLOPSPERCALC 2

float fa[FLOPS_ARRAY_SIZE] __attribute__((aligned(64)));
float fb[FLOPS_ARRAY_SIZE] __attribute__((aligned(64)));

int main(int argc, char *argv[]) {
	int i, j, k;
	int numthreads;
	double tstart, tstop, ttime;
	double gflops = 0.0;
	float a = 1.1;

	printf("Initialising.\n");

	omp_set_num_threads(122);
	kmp_set_defaults("KMP_AFFINITY=scatter");
	
	#pragma omp parallel for
	for (i = 0; i < FLOPS_ARRAY_SIZE; i++) {
		if (i == 0) {
			numthreads = omp_get_num_threads();
		}
		fa[i] = (float) i + 0.1;
		fb[i] = (float) i + 0.2;
	}

	printf("Starting Compute for %d Threads.\n", numthreads);
	
	tstart = dtime();

	#pragma omp parallel for private(j,k)
	for (i = 0; i < numthreads; i++) {
		int offset = i * LOOP_COUNT;
		for (j = 0; j < MAXFLOPS_ITERS; j++) {
			for (k = 0; k < LOOP_COUNT; k++) {
				fa[k + offset] = a * fa[k + offset] + fb[k + offset];
			}
		}
	}
	
	tstop = dtime();

	gflops = (double) (1.0e-9 * LOOP_COUNT * MAXFLOPS_ITERS * FLOPSPERCALC);

	ttime = tstop - tstart;

	if (ttime > 0.0) {
		printf("GFlops = %10.3lf, Secs = %10.3lf, GFlops/sec = %10.3lf\r\n",
				gflops, ttime, gflops/ttime);
	}

	return 0;
}
