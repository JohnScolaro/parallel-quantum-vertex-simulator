set term gif
set xrange [-4:4]
set yrange [-2:3.5]
set key below
set style data line
set title "#YOLO"
set xlabel "Position"
set ylabel "Probability"

plot "data" using 1:2 title 'Probability', "data" using 1:3 title 'Real', "data" using 1:4 title 'Imaginary'
