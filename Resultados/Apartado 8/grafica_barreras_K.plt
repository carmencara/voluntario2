fichero = "barreras_K.txt"

salida = "barreras_K.png"

set xlabel "1/V (V^{-1})" 
set ylabel "1/r^2 (m^{-2})"


f(x)=a*exp(-b*x)
fit f(x) fichero via a,b

#set yrange [-1.8:0]
#set ytics (-1.8,-1.6,-1.4,-1.2,"-1.0"-1,-0.8,-0.6,-0.4,-0.2,0)

unset key

plot fichero w xyerrorbars lc rgb "red" ps 0.2 pt 7

set style lines 1 lt rgb "pink" 
replot f(x) with lines linestyle 1 dt 1 lw 2


set term pngcairo enhanced color font "Arial,14"
set output salida
replot
