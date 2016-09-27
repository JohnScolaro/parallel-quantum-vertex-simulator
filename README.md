#Introduction

This is a project to create some well documented code to evolve a 1 dimensional wave-function according to the Schroedinger Equation via a 'Forwards in Time, Central Difference' (FCTS) approach to the Crank-Nicolson method. It utilizes OpenMP to simplify multiprocessing for the Intel Xeon Phi, which should lead to increased computation speeds.

###Overview of files in this folder:

*README.md:*
This file.

*quantum_sim.c:*
The file which contains all the computation.

*plotter.gnuplot:*
This is a script used to custimize the style of gnuplot which is used to plot our data, and later create animations.

*makemove.sh:*
This script calls gnuplot and convert to take out series of data outputs, and turn them into a series of gifs, and finally an animation.

*Makefile:*
This makes the quantum_sim.

*plots:*
This folder will contain all the data from quantum_sim when it is a ran. If it doesn't exist, quantum_sim will refuse to run until one is created in the same folder as it.

#Instructions
Simply type: 'Make'. This will make the quantum_sim program. Then run it with './quantum_sim'. This generated all the data in the plots folder. Then run './make_movie.sh'. This will create a gif in the plots folder, which will be an animation of the wavefunction over time. The get the program to create different plots, or to do custom things to the animation, you can either edit the variable values at the top of quantum_sim.c, and recompile and run the code again, or you can add your own code to quantum_sim.c, and rerun it.
