#!/usr/bin/gnuplot
reset
set terminal png

set xrange [1:4]
set xlabel "qgram size"

set ylabel "time (seconds)"
#set yrange [100.0:200.0]
set yrange [50:250]

set title "Inverted - DNA preprocessing"
#set key reverse Left outside
set grid

set style data linespoints

plot "inverted/chr1/summary.plot" using 1:2 title "c1, e5", \
"" using 1:3 title "c2, e5"
#
