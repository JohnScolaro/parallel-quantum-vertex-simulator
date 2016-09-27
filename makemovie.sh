#!/bin/bash
cd plots/;
for i in dat*.dat;
	do cp $i data;
	gnuplot "../plotter.gnuplot" > $i.gif 2> /dev/null;
	done
rm data
convert -delay 5 dat*.gif animation.gif
#rm dat*.gif
cd ../;
