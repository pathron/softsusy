set pm3d map
set terminal png
set xlabel "m0(GeV)"
set ylabel "m12(GeV)"
set palette defined (0 "purple", 12.5 "blue", 25 "#00ffff", 37.5 "green", 50 "yellow", 62.5 "orange", 75 "red", 87.5 "#990000", 100 "black", 101 "grey")
set output "DelBG.png"
set title "Del_BG"
set zlabel "Del-BG"
splot "priorscan.dat" using 1:2:16 with points pointtype 7 pointsize 3 palette
set zlabel "Del-BG-noYt"
set title "Del-BG-noYt"
set output "DelBG_noYt.png"
splot "priorscan.dat" using 1:2:17 with points pointtype 7 pointsize 3 palette


