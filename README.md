#Introduction

This is a project to create some well documented code to evolve a 1 dimensional wave-function according to the Schroedinger Equation via a 'Forwards in Time, Central Difference' (FCTS) approach to the Crank-Nicolson method. It utilizes OpenMP to simplify multiprocessing for the Intel Xeon Phi, which should lead to increased computation speeds.

###Overview of files in this folder:

####README.md:
This file.

####quantum_sim.c:
The file which contains all the computation.

####plotter.gnuplot:
This is a script used to custimize the style of gnuplot which is used to plot our data, and later create animations.

####makemove.sh:
This script calls gnuplot and convert to take out series of data outputs, and turn them into a series of gifs, and finally an animation.

####Makefile:
This makes the quantum_sim.

####plots:
This folder will contain all the data from quantum_sim when it is a ran. If it doesn't exist, quantum_sim will refuse to run until one is created in the same folder as it.

####Examples:
This folder contains some example gifs for viewing. It hosts the images which displayed in this README file.

#Instructions
Simply type: 'Make'. This will make the quantum_sim program. Then run it with './quantum_sim'. This generates all the data in the plots folder. Then run './make_movie.sh'. This will create a gif in the plots folder, which will be an animation of the wavefunction over time. The get the program to create different plots, or to do custom things to the animation, you can either edit the variable values at the top of quantum_sim.c, and recompile and run the code again, or you can add your own custom code to quantum_sim.c, and rerun it.

#Example Options
Some quick examples to get you going immediately:
###1: The Ground State of a Harmonic Oscillator
	g = 0.0;
	width = 1;
	choice = 0;
	xinit = 0.0;
	k = 0.0;
	x0 = -6.0;
	xmax = 6.0;
	Nx = 600;
	dt = 0.05;
	Nt = 100;
	normalize = 1;
	moving = 0;
	moveSpeed = 0.5;
![](https://github.com/JohnScolaro/parallel-quantum-vertex-simulator/blob/master/examples/ground_state.gif)
###2: A displaced third Eigenstate of the Harmonic Oscillator evolving through Negative Imaginary Time
	g = 0.0;
	width = 1;
	choice = 3;
	xinit = 0.6;
	k = 0.0;
	x0 = -6.0;
	xmax = 6.0;
	Nx = 600;
	dt = -0.03 * I;
	Nt = 200;
	normalize = 1;
	moving = 0;
	moveSpeed = 0.5;
![](https://github.com/JohnScolaro/parallel-quantum-vertex-simulator/blob/master/examples/neg_imag.gif)
###3: The Second Lowest Eigenstate, with a Forced Moving Zero and Smoothing
	g = 0.0;
	width = 1;
	choice = 2;
	xinit = 0.0;
	k = 0.0;
	x0 = -6.0;
	xmax = 6.0;
	Nx = 600;
	dt = 0.05;
	Nt = 200;
	normalize = 1;
	moving = 1;
	moveSpeed = 0.3;
![](https://github.com/JohnScolaro/parallel-quantum-vertex-simulator/blob/master/examples/forced_zero.gif)
