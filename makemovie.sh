#!/bin/bash
cd plots/;
for i in dat*.dat;
	do cp $i data;
	gnuplot "../plotter.gnuplot" > $i.gif;
	done
rm data
convert dat*.dat animation.gif
rm dat*.gif
cd ../;
