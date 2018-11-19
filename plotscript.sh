#!/bin/bash

SIM="detonation"
FOLDERPATH="./output/$SIM"

mkdir -p "$FOLDERPATH/plots"


gnuplot<<EOF
reset
unset key
set xtics 0.5
set terminal pdfcairo size 20,10 font 'times,16' rounded

set title "$SIM: p-v space"
set xlabel "specific volume (v)"
set ylabel "pressure (p)"
set xrange [0:4]
set style line 1 lc 7 lw 3 pt 7 ps 0.5

set arrow from 1.12,1.05 to 0.95,1.5 lw 6
set label "time" at 0.83,1.4

filenums="102 105 109 110 111"

set output "$FOLDERPATH/plots/Rlines1.pdf"
plot for [n in filenums] "$FOLDERPATH/snapshot_".n.".dat" u 5:4 w p ls 1

unset arrow
unset label

set output "$FOLDERPATH/plots/Rlines2.pdf"
plot "$FOLDERPATH/snapshot_118.dat" u 5:4 w p ls 1

set xlabel 'position'
set xtics 0.1
unset ylabel
set xrange [0:1]
#set arrow from 0.06,3.4 to 0.23,3.4
#set label "time" at 0.03,3.42

filenums="75 90 105 106 107 108 109 110 111 112 115 122 129"

set title '$SIM: Density'
set output "$FOLDERPATH/plots/density.pdf"
plot for [n in filenums] "$FOLDERPATH/snapshot_".n.".dat" u 1:2 w l ls 1

set title '$SIM: Pressure'
set output "$FOLDERPATH/plots/pressure.pdf"
plot for [n in filenums] "$FOLDERPATH/snapshot_".n.".dat" u 1:4 w l ls 1

set title '$SIM: Temperature' font 'times,12'
set output "./output/$SIM/plots/temperature.pdf"
plot for [n in filenums] "$FOLDERPATH/snapshot_".n.".dat" u 1:7 w l ls 1

set title '$SIM: Species' font 'times,12'
set output "./output/$SIM/plots/species.pdf"
plot for [n in filenums] "$FOLDERPATH/snapshot_".n.".dat" u 1:8 w l ls 1
EOF
