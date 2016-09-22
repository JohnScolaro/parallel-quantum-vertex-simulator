set term gif
set xrange [-3:3]
set yrange [0:2]
set key below
set style data line
set title "G < 0"
set xlabel "Position"
set ylabel "Probability"

plot "data"
