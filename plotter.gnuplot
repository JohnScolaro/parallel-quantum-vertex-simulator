set term gif
set xrange [-1:1]
set yrange [0:2]
set key below
set style data line
set title "Progress Report"
set xlabel "Position"
set ylabel "Probability"

plot "data"
