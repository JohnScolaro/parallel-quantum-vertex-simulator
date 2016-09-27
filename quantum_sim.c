/*
 * This program solves the Schroedinger Equations by using a forwards in time,
 * central difference Crank-Nicolson solver. It has the ability to use OpenMP
 * to easily parallelize the solution, and solve the equation much faster with
 * a larger amount of computational cores.
 */


/* Includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <complex.h>


/* Defines */
#define PI 3.1415926


/* Function Predeclarations */
int tridag(double complex *, double complex *, double complex *, double complex *, double complex *, int);
double func(double, double, double, double, int);


int main (int argc, char *argv[]) {
	// norm -> Norm of the solution
	// avloc -> average location
	// name -> name of the file containing norm and expectation values over time
	// probdens -> probability density
	// reamp -> real part of the amplitude
	// imamp -> imaginary part of the amplitude
	// alpha ->	collection of constants in the schroedinger equation, assuming
	//			that hbar, weight, and mass are all equal to 1.
	// beta -> collection of other factors
	// var -> variance, and var is used to calculate this


	/* Declarations */
	double width, x0, xmax, x, xinit, dx, t = 0.0, dt;
	double norm = 0.0, avloc = 0.0;
	double var = 0.0, var1 = 0.0;
	double probdens, reamp, imamp, g, k;
	double complex alpha, beta;
	int Nx, Nt, n = 0, j, ierr, choice, moving;


	/*
	 * Constant setting:
	 * If you want to run the program with different initial parameters, change
	 * these numbers.
	 */

	/* Name of output file */
 	char name[20] = "test.txt";
	/* Non-linear coupling constant */
	g = 0.0;
	/* Width of initial function */
	width = 1.0;
	/* Initial function. 0 = Gaussian, 1 = hyperbolic secant */
	choice = 2;
	/* The position of the initial function */
	xinit = 0.0;
	/* The momentum kick of the initial waveform */
	k = 0.0;
	/* Position of the first spacial grid boundary (x_min) */
	x0 = -6.0;
	/* Position of the second spacial grid boundary (x_max) */
	xmax = 6.0;
	/* The number of spacial grid points */
	Nx = 600;
	/* Time Step */
	dt = 0.05;
	/* Number of time steps */
	Nt = 400;
	/* Option to force a moving zero in the wavefunction. 0 = no, 1 = yes. */
	moving = 1;

	/*
	 * Calculated Variables
	 */
	/* The distance between two points. */
	dx = (xmax - x0) / ((double) Nx + 1.0);

	beta = 0.5 * dt * I;
	alpha = -0.5 * beta / (dx * dx);

	/*
	 * Allocating space for arrays to store values in to compute the tridiagonal
	 * matrix.
	 */
	double complex *matxold, *adiag, *alower, *aupper, *rrhs, *xsoln, *potxtold;
	double complex *potxtnew;
	matxold = malloc(Nx * sizeof(double complex));
	adiag = malloc(Nx * sizeof(double complex));
	alower = malloc(Nx * sizeof(double complex));
	aupper = malloc(Nx * sizeof(double complex));
	rrhs = malloc(Nx * sizeof(double complex));
	xsoln = malloc(Nx * sizeof(double complex));
	potxtold = malloc(Nx * sizeof(double complex));
	potxtnew = malloc(Nx * sizeof(double complex));

	/* Set up files */
	FILE *fp;
	FILE *fp2;
	fp2 = fopen(name, "w");
	fprintf(fp2, "%15s%15s%15s%20s%15s\n", "Iteration", "Time (s)", "Norm", "Av. location (m)", "Variance");

	/*
	 * Calculate the initial values of the wave function. This is done by
	 * setting matxold values to the output of func for a given 'choice'
	 * variable.
	 */
	x = x0;
	for (j = 0; j < Nx; j++) {
		x = x + dx;
		matxold[j] = func(width, x, xinit, k, choice);
		probdens = cabs(matxold[j]) * cabs(matxold[j]);
		potxtold[j] = 0.5 * pow(x, 2.0); // g * probdens once I get it working for an oscillator
		norm = norm + dx * probdens; //The norm
		avloc = avloc + dx * x * probdens; // The average position
		var1 = var1 + dx * x * x * probdens; // Not super important
	}
	var = pow((var1 - pow(avloc, 2.0)), 0.5);
	fprintf(fp2, "%15.5i %15.5lf %15.5lf %15.5lf %15.5lf\n", n, t, norm, avloc, var);
	x = x0;

	/* Put the initial values into the first file. */
	fp = fopen("plots/dat0000.dat", "w");
	if (fp == NULL) {
		printf("This program won't create data unless you have a plots"
				"directory in the same directory as this file. To fix, type:"
				" \"mkdir plots\".\n");
		exit(-1);
	}

	for (j = 0; j < Nx; j++) {
		x = x + dx;
		probdens = cabs(matxold[j]) * cabs(matxold[j]);
		fprintf(fp, "%lf %lf %lf %lf\n", x, probdens, creal(matxold[j]), cimag(matxold[j]));
	}
	fclose(fp);

	/* For every timestep */
	for (n = 1; n <= Nt; n++) {
		t = t + dt;
		norm = 0.0;
		avloc = 0.0;
		var1 = 0.0;

		/* Generate a file with a unique name */
		char fname[64];
		strcpy(fname, "plots/datxxxx");
		sprintf((fname + 9), "%04i", n);
		strcat(fname, ".dat");
		fp = fopen(fname, "w");

		/*
		 * For every point in space, build up the A and R.H.S matricies for
		 * the tridiagonal solver.
		 */
		for (j = 0; j < Nx; j++) {
			if (j == 0) {
				rrhs[j] = (1.0 + (2.0 * alpha) - (beta * potxtold[j])) *
						matxold[j] - (alpha * matxold[j + 1]);
			} else if (j == Nx - 1) {
				rrhs[j] = -alpha * matxold[j - 1] + (1.0 + 2.0 * alpha - beta *
						potxtold[j]) * matxold[j];
			} else {
				rrhs[j] = -alpha * matxold[j - 1] + (1.0 + 2.0 * alpha - beta *
						potxtold[j]) * matxold[j] - alpha * matxold[j + 1];
			}
			probdens = cabs(matxold[j]) * cabs(matxold[j]);
			potxtnew[j] = potxtold[j] - g * probdens;
			alower[j] = alpha;
			adiag[j] = (1.0 - 2.0 * alpha + beta * potxtnew[j]);
			aupper[j] = alpha;
		}

		/* Send to tridiagonal solver */
		ierr = tridag(adiag, alower, aupper, rrhs, xsoln, Nx);

		/* Check for errors */
		if (ierr == -1) {
			printf("ierr = -1 so no tridag solution found!\n");
			exit(-1);
		}

		/* Manually edit the waveform to attempt to make it travel. Edit all
		 * points on either side to smooth the waveform. The point of movement
		 * is: 300 + (n / 4), and it is smoothed exponentially on either side.
		 */
		if (moving == 1) {
			for (int j = 0; j < Nx; j++) {
				if (j == 300 + (n / 4)) {
					// Fixing the centre point
					xsoln[j] = 0.0;
				} else if (j < (300 + (n / 4))) {
					// Smoothing below the centre point
					xsoln[j] = (1 - exp(-0.2 * ((300 + (n / 4)) - j))) * xsoln[j];
				} else if (j > (300 + (n / 4))) {
					// Smoothing above the centre point
					xsoln[j] = (1 - exp(-0.2 * (j - (300 + (n / 4))))) * xsoln[j];
				}
			}
		}

		/*
		 * Run through the solution array to calculate monitored values at the
		 * current time and print them.
		 */
		x = x0;
		for (j = 0; j < Nx; j++) {
			x = x + dx;
			probdens = cabs(xsoln[j]) * cabs(xsoln[j]);
			reamp = creal(xsoln[j]);
			imamp = cimag(xsoln[j]);
			fprintf(fp, "%lf %lf %lf %lf\n", x, probdens, reamp, imamp);
			norm = norm + dx * probdens;
			avloc = avloc + dx * x * probdens;
			var1 = var1 + dx * x * x * probdens;
		}
		var = pow(var1 - pow(avloc, 2.0), 0.5);
		fprintf(fp2, "%15.5i %15.5lf %15.5lf %15.5lf %15.5lf\n", n + 1, t, norm,
		 		avloc, var);
		memcpy(matxold, xsoln, Nx * sizeof(double complex));
		fclose(fp);
	}
	fclose(fp2);

	/* Free allocated memory */
	free(matxold);
	free(adiag);
	free(alower);
	free(aupper);
	free(rrhs);
	free(xsoln);
	free(potxtold);
	free(potxtnew);

	return 0;
}


/* Tridiagonal solver from Michael */
int tridag(double complex adiag[], double complex alower[], double complex aupper[], double complex rrhs[], double complex xsoln[], int n) {
	int j;
	double complex bet, *gam;

	/* Error checking */
	if (fabs(adiag[0] == 0.0)) {
		return -1;
	}

	gam = malloc(n * sizeof(double complex));

	bet = adiag[0];
	xsoln[0] = rrhs[0] / bet;

	for (j = 1; j < n; j++){
		gam[j] = aupper[j - 1] / bet;
		bet = adiag[j] - alower[j] * gam[j];
		if (fabs(bet) == 0.0) {
			return(-1);
		}
		xsoln[j] = (rrhs[j] - alower[j] * xsoln[j - 1]) / bet;
	}

	for (j = n - 2; j >= 0; j--) {
		xsoln[j] -= gam[j + 1] * xsoln[j + 1];
	}

	free(gam);

	return 0;
}


/*
 * This function creates initial functions depending on what the value in
 * 'choice' is.
 *
 * 0 = Gaussian
 * 1 = Hyperbolic Secant
 * 2 = Second Energy Level Eigenstate of a Harmonic Oscillator
 * 3 = Third Energy Level Eigenstate of a Harmonic Oscillator
 */
double func(double width, double x, double xinit, double k, int choice) {
	double complex ans;

	/* Gaussian */
	if (choice == 0) {
		ans = (1.0 / pow(width * pow(PI, 0.5), 0.5)) * exp(-pow(x - xinit, 2.0)
				/ (2.0 * width * width));
	/* Hyperbolic Secant */
	} else if (choice == 1) {
		ans = pow(PI / (6.0 * width * pow(2.0, 0.5)), 0.5) *
				(1.0 / cosh((PI / (3.0 * pow(2.0, 0.5))) *
				((x - xinit) / width))) * cexp(I * k * x);
	/* Second Eigenstate of a Harmonic Oscillator */
	} else if (choice == 2) {
		ans = (1.0 / pow(width * pow(PI, 0.5), 0.5)) * exp(-pow(x - xinit, 2.0)
				/ (2.0 * width * width)) * 2 * (x - xinit);
	} else if (choice == 3) {
			ans = (1.0 / pow(width * pow(PI, 0.5), 0.5)) * exp(-pow(x - xinit, 2.0)
					/ (2.0 * width * width)) * (4.0 * pow(x - xinit, 2.0) - 2.0);
	}

	return ans;
}
