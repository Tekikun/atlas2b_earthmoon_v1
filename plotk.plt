set xlabel "a (ER)"
set ylabel "streng"
set logscale y 10
set yrange [1e-5:0.02]
set xrange [5:80]
plot "atlas2bres.dat" u 4:9 with impulses lw 1
pause -1

