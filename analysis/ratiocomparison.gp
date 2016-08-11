plot 'no-sample-1801.txt' using 1:2 w l lw 3 title 'Corrected with no-sample one-arm scan and sqrt'
replot 'sample-1801.txt' using 1:2 w l lw 3 title 'Corrected with sample one-arm scan and sqrt'
replot 'no-sqrt-1801.txt' using 1:2 w l lw 3 title 'Corrected with sample one-arm scan and no sqrt'
replot 'no-sample-sqrt-1801.txt' using 1:2 w l lw 3 title 'Corrected with no-sample one-arm scan and no sqrt'
set grid

set terminal postscript color font 'courier' 
set output '18-01-2 150 GHz Xpol ratio comparison.ps'
replot
