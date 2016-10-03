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
#include <omp.h>


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
	// reamp -> real part of psi
	// imamp -> imaginary part of psi
	// alpha ->	collection of constants in the schroedinger equation, assuming
	//			that hbar, weight, and mass are all equal to 1.
	// beta -> same as alpha
	// var -> variance, and var is used to calculate this


	/* Declarations */
	double width, x0, xmax, x, xinit, dx;
	double norm = 0.0, avloc = 0.0;
	double var = 0.0, var1 = 0.0;
	double probdens, reamp, imamp, g, k, moveSpeed;
	double complex alpha, beta, dt, t = 0.0;
	int Nx, Nt, n = 0, j, ierr, choice, moving, output_number = 1;
	int output_frequency, normalize;


	/*
	 * Constant setting:
	 * If you want to run the program with different initial parameters, change
	 * these numbers.
	 */

	/* Name of output file */
 	char name[20] = "test";
	/* Non-linear coupling constant */
	g = 0.0;
	/* Width of initial function */
	width = 1;
	/* Initial function. Check the 'func' function to see all options.*/
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
	/* Time Step. Simply add (-1 * I) to make it imaginary */
	dt = 0.05;
	/* Number of time steps */
	Nt = 200;
	/* Option to normalize the probability after each time step */
	normalize = 1;
	/* Option to force a moving zero in the wavefunction. 0 = no, 1 = yes. */
	moving = 1;
	/* Option to change the moving speed */
	moveSpeed = 0.3;
	/*
	 * How often the program prints data. 1 = every time step, 2 = every 2nd,
	 * etc. This is useful for longer, more precise calculations, as printing
	 * requires way more overhead than simply performing arithmatic, so for
	 * very precise solutions it may be useful to only print every 50 or 100
	 * time steps.
	 */
	output_frequency = 1;

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
	double complex *psi, *adiag, *alower, *aupper, *rrhs, *xsoln, *potential;
	double complex *potentialInteraction;
	/* This is psi */
	psi = malloc(Nx * sizeof(double complex));
	/* This is the potential */
	potential = malloc(Nx * sizeof(double complex));
	potentialInteraction = malloc(Nx * sizeof(double complex));
	/* These are arrays populated later, and used in the tridiagonal solver */
	adiag = malloc(Nx * sizeof(double complex));
	alower = malloc(Nx * sizeof(double complex));
	aupper = malloc(Nx * sizeof(double complex));
	rrhs = malloc(Nx * sizeof(double complex));
	xsoln = malloc(Nx * sizeof(double complex));

	/* Set up files */
	FILE *fp = NULL, *fp2 = NULL;
	fp2 = fopen(name, "w");
	fprintf(fp2, "%15s%15s%15s%20s%15s\n", "Iteration", "Time", "Norm",
			"Av. location", "Variance");

	/*
	 * Calculate the initial values of the wave function. This is done by
	 * setting psi values to the output of func for a given 'choice'
	 * variable.
	 */
	x = x0;
	for (j = 0; j < Nx; j++) {
		x = x + dx;
		psi[j] = func(width, x, xinit, k, choice);
		probdens = cabs(psi[j]) * cabs(psi[j]);
		potential[j] = 0.5 * pow(x, 2.0); // This is the potential for a harmonic oscillator
		norm += dx * probdens; //The norm
		avloc += dx * x * probdens; // The average position
		var1 += dx * x * x * probdens; // Not super important
	}
	var = pow((var1 - pow(avloc, 2.0)), 0.5);
	fprintf(fp2, "%15.5i %15.5lf %15.5lf %15.5lf %15.5lf\n", n, t, norm, avloc, var);
	x = x0;

	/* Put the initial values into the first file. */
	fp = fopen("plots/dat0000.dat", "w");
	if (fp == NULL) {
		printf("This program won't create data unless you have a plots "
				"directory in the same directory as this file. To fix, type: "
				"\"mkdir plots\".\n");
		exit(-1);
	}

	for (j = 0; j < Nx; j++) {
		x = x + dx;
		probdens = cabs(psi[j]) * cabs(psi[j]);
		fprintf(fp, "%lf %lf %lf %lf\n", x, probdens, creal(psi[j]),
				cimag(psi[j]));
	}
	fclose(fp);

	/*
	 * This is where the program essentially begins doing the bulk of its work
	 * This first loop loops over each timestep.
	 */
	for (n = 1; n <= Nt; n++) {
		t = t + dt;

		/* Generate a file with a unique name. */
		if (n % output_frequency == 0) {
			char fname[64];
			strcpy(fname, "plots/datxxxx");
			sprintf((fname + 9), "%04i", output_number);
			strcat(fname, ".dat");
			fp = fopen(fname, "w");
			output_number++;
		} else {
			fp = NULL;
		}


		/*
		 * For every point in space, build up the A and R.H.S matricies for
		 * the tridiagonal solver. In this for loop, we update rrhs depending
		 * on the values of psi (psi), constants, and the potential.
		 * Because of this, we can parallelize it with OpenMP.
		 */
		//omp_set_num_threads(20);
		#pragma omp parallel for
		for (j = 0; j < Nx; j++) {
			if (j == 0) {
				rrhs[j] = (1.0 + (2.0 * alpha) - (beta * potential[j])) *
						psi[j] - (alpha * psi[j + 1]);
			} else if (j == Nx - 1) {
				rrhs[j] = -alpha * psi[j - 1] + (1.0 + 2.0 * alpha - beta *
						potential[j]) * psi[j];
			} else {
				rrhs[j] = -alpha * psi[j - 1] + (1.0 + 2.0 * alpha - beta *
						potential[j]) * psi[j] - alpha * psi[j + 1];
			}
			probdens = cabs(psi[j]) * cabs(psi[j]);
			potentialInteraction[j] = potential[j] - g * probdens;
			alower[j] = alpha;
			adiag[j] = (1.0 - 2.0 * alpha + beta * potentialInteraction[j]);
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
		 * is: ((Nx / 2) + (int) ((Nx * n * moveSpeed) / (2 * Nt))), and it is
		 * smoothed exponentially on either side.
		 *
		 * This formula lets us start in the middle, at Nx / 2, and snowly
		 * travel to (Nx / 2) + (Nx * moveSpeed / 2) by the end of the sim. So
		 * if moveSpeed equals 1, we will move from Nx / 2 to Nx. If Movespeed
		 * -0.5, then it will move from Nx / 2 to Nx / 4 over the whole sim.
		 */
		if (moving == 1) {
			#pragma omp parallel for
			for (int j = 0; j < Nx; j++) {
				if (j == ((Nx / 2) + (int) ((Nx * n * moveSpeed) / (2 * Nt)))) {
					// Fixing the centre point
					xsoln[j] = 0.0;
				} else if (j <= ((Nx / 2) + (int) ((Nx * n * moveSpeed) /
						(2 * Nt)))) {
					// Smoothing below the centre point
					xsoln[j] *= (1 - exp(-0.2 * (((Nx / 2) + (int)
							((Nx * n * moveSpeed) / (2 * Nt))) - j)));
				} else if (j >= ((Nx / 2) + (int) ((Nx * n * moveSpeed) /
						(2 * Nt)))) {
					// Smoothing above the centre point
					xsoln[j] *= (1 - exp(-0.2 * (j - ((Nx / 2) + (int)
							((Nx * n * moveSpeed) / (2 * Nt))))));
				}
			}
		}

		/*
		 * Run through the solution and calculate monitored values at the
		 * current time and print them. Also prints the values which we graph.
		 */
		x = x0;
		norm = 0.0;
		avloc = 0.0;
		var1 = 0.0;
		for (j = 0; j < Nx; j++) {
			x = x + dx;
			probdens = cabs(xsoln[j]) * cabs(xsoln[j]);
			reamp = creal(xsoln[j]);
			imamp = cimag(xsoln[j]);
			if (fp != NULL) {
				fprintf(fp, "%lf %lf %lf %lf\n", x, probdens, reamp, imamp);
			}
			norm += dx * probdens;
			avloc += dx * x * probdens;
			var1 += dx * x * x * probdens;
		}

		var = pow(var1 - pow(avloc, 2.0), 0.5);
		fprintf(fp2, "%15.5d %15.5lf %15.5lf %15.5lf %15.5lf\n", n, t, norm,
		 		avloc, var);

		/*
		 * Normalize the wavefunction after each time step.
		 */
		if (normalize == 1) {
			#pragma omp parallel for
			for (j = 0; j < Nx; j++) {
				xsoln[j] /= sqrt(norm);
			}
		}

		memcpy(psi, xsoln, Nx * sizeof(double complex));
		if (fp != NULL) {
			fclose(fp);
		}
	}
	fclose(fp2);

	/* Free allocated memory */
	free(psi);
	free(adiag);
	free(alower);
	free(aupper);
	free(rrhs);
	free(xsoln);
	free(potential);

	return 0;
}


/* Tridiagonal solver from Michael Bromley*/
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

	for (j = 1; j < n; j++) {
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
