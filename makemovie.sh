#!/bin/bash
cd plots/;
for i in dat*.dat;
	do cp $i data;
	gnuplot "../plotter.gnuplot" > $i.gif;
	done
rm data
convert -delay 2 dat*.gif animation.gif
rm dat*.gif
cd ../;
