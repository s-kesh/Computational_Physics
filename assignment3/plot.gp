#set term pngcairo size 1366, 768 enhanced
#set out "explicit_1d.png"
#set title "Explicit 1D"
#set title "Implicit 1D"
set title "Crank Nicolson 1D"
set xrange [0 : 1]
set yrange [0 : 15]
set xlabel "Space Interval"
#set ylabel "Time step (\del t = 0.0033)"
set ylabel "Time step (\del t = 0.05)"
#set cbrange [0 : 1]

set table 'test.dat'
splot "data4" matrix u ( 0.05 * $1):2:3 notitle
unset table

set contour base
set cntrparam level increment 0, 0.1, 1
unset surface
set table 'cont.dat'
splot "data4" matrix u ( 0.05 * $1):2:3 notitle
unset table

reset
set xrange [0 : 1]
set yrange [0 : 15]
set xlabel "Space Interval"
#set ylabel "Time step (\del t = 0.0033)"
set ylabel "Time step (\del t = 0.05)"
unset key
set palette rgbformulae 22,13,10
plot 'test.dat' with image notitle, 'cont.dat' w l lt -1 lw 0.5 notitle
