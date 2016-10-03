set term gif
set xrange [-3:3]
set yrange [0:1]
set key below
set style data line
set title "Evolving through Negative Imaginary Time"
set xlabel "Position"
set ylabel "Normalized Probability"

#Used to plot the probability, the real and imaginary axes of psi
#plot "data" using 1:2 title 'Probability', "data" using 1:3 title 'Real', "data" using 1:4 title 'Imaginary'

#Used to plot the probability only
plot "data" using 1:2 title 'Probability'
