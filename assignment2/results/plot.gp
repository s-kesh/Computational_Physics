set term pdf
set out "Method_Comparision.pdf"
set xlabel "t"
set ylabel "y"
set title "Method Comparision (step size = 0.1)"
plot "eular_0.1.non.nn" w lp t "Euler Method",  "avg_0.1.non.nn" w lp t "Average Method", "mid_0.1.non.nn" w lp t "MidPoint Method", "rk4_0.1.non.nn" w lp t "Rounge Kutta 4 Method", exp(-0.5 * x**2) w l t "analytical solution"
