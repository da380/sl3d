reset
unset key

set pm3d interpolate 2,2
set view map

unset xtics
unset ytics

#set cbrange [-1200:200]


# COLORBOX
#set terminal postscript eps enhanced color font "Times-Roman,20" size 6,0.75
#set colorbox horizontal user origin 0.02,0.7 size 0.96, 0.2
#set output 'sphvisccb2.eps'

#set cblabel "Change in sea level/ m"

#set lmargin at screen -1.0
#set rmargin at screen -1.0
#set tmargin at screen 1.5
#set bmargin at screen 1.5


# PLOT
set terminal pngcairo enhanced font "Times-Roman,20" size 1200,600

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]

unset colorbox

set lmargin at screen 0.02
set rmargin at screen 0.99
set tmargin at screen 0.99
set bmargin at screen 0.02

#ACTUAL SL - 1D
set output 'sl_test.png'




splot 'sl.64.0.0' u 4:5:7 w pm3d, 'coast_lines_110_rob.dat' u 1:2:(0) w l lw 3 lc rgb "black" notitle



