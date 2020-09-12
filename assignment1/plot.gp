set term pdf size 5, 12
set key outside
set size ratio 3
set out "plot_all.pdf"
set sample 1000000
set grid

# Potential Well
f(x) = x > -10 && x < 10 ? 0 : 10

# Eigen Functions
f0(x) = x < -10 ? 27934512487.75554367192456266980 * exp(16113944159.97361865523000000000 * 1e-10 * x) :	(x < 10 ? 30684.98192921710832648474 * cos(1479253369.84386450968000000000 * 1e-10 * x) : 27934512487.75554367192456266980 * exp(-16113944159.97361865523000000000 * 1e-10 * x))
f1(x) = x < -10 ? -45491522587.91441059191056051883 * exp(15909084146.12287410487000000000 * 1e-10 * x) :	(x < 10 ? 30673.44440065669395103787 * sin(2957774256.85668941029000000000 * 1e-10 * x) : -1 * -45491522587.91441059191056051883 * exp(-15909084146.12287410487000000000 * 1e-10 * x))
f2(x) = x < -10 ? -48180845119.17715227526551124513 * exp(15562138353.74040447933000000000 * 1e-10 * x) :	(x < 10 ? 30653.24327255077877820927 * cos(4434775843.48794488980000000000 * 1e-10 * x) : -48180845119.17715227526551124513 * exp(-15562138353.74040447933000000000 * 1e-10 * x))
f3(x) = x < -10 ? 38977320305.05660330132450563282 * exp(15064094751.96157081591000000000 * 1e-10 * x) :	(x < 10 ? 30622.69372807396857598458 * sin(5909351591.06870030584000000000 * 1e-10 * x) : -1 * 38977320305.05660330132450563282 * exp(-15064094751.96157081591000000000 * 1e-10 * x))
f4(x) = x < -10 ? 25036767307.21036535596817554425 * exp(14400603105.44173050420000000000 * 1e-10 * x) :	(x < 10 ? 30578.87306617481453971213 * cos(7380380554.04220380509000000000 * 1e-10 * x) : 25036767307.21036535596817554425 * exp(-14400603105.44173050420000000000 * 1e-10 * x))
f5(x) = x < -10 ? -12786580105.53423152434312234550 * exp(13549505440.06645974526000000000 * 1e-10 * x) :	(x < 10 ? 30516.70266388035250220034 * sin(8846371530.32539848049000000000 * 1e-10 * x) : -1 * -12786580105.53423152434312234550 * exp(-13549505440.06645974526000000000 * 1e-10 * x))
f6(x) = x < -10 ? -5076313223.74509295086111310901 * exp(12476003773.01723361417000000000 * 1e-10 * x) :	(x < 10 ? 30426.86354082132231935677 * cos(10305179124.04253508747000000000 * 1e-10 * x) : -5076313223.74509295086111310901 * exp(-12476003773.01723361417000000000 * 1e-10 * x))
f7(x) = x < -10 ? 1488643900.57666701414840364482 * exp(11122284287.27891039155000000000 * 1e-10 * x) :	(x < 10 ? 30290.38290918199107156658 * sin(11753390113.32011304929000000000 * 1e-10 * x) : -1 * 1488643900.57666701414840364482 * exp(-11122284287.27891039155000000000 * 1e-10 * x))
f8(x) = x < -10 ? 290602828.26898421280301370865 * exp(9381296217.93909780004000000000 * 1e-10 * x) :	(x < 10 ? 30061.15292712211726786292 * cos(13184789273.78909227290000000000 * 1e-10 * x) : 290602828.26898421280301370865 * exp(-9381296217.93909780004000000000 * 1e-10 * x))
f9(x) = x < -10 ? -29495396.59862413989205762755 * exp(7008679092.53200717965000000000 * 1e-10 * x) :	(x < 10 ? 29582.68858399947625951665 * sin(14585122704.34661238368000000000 * 1e-10 * x) : -1 * -29495396.59862413989205762755 * exp(-7008679092.53200717965000000000 * 1e-10 * x))
f10(x) = x < -10 ? -555588.15380842822015622367 * exp(3026693282.95552197975000000000 * 1e-10 * x) :	(x < 10 ? 27416.36866337119577694642 * cos(15896116339.96978030826000000000 * 1e-10 * x) : -555588.15380842822015622367 * exp(-3026693282.95552197975000000000 * 1e-10 * x))

# Makeup
set title "Eigen Functions of Bound States.\nY-axis of Eigenfunctions have been scaled down to fit in plot.\nSo, do not scale."
set xrange [-18 : 18 ]
set yrange [-2 : 13]
set xlabel "Width (A⁰)"
set ylabel "Potential (eV)"
set ytics (0.08357, 0.33104, 0.75110, 1.33362, 2.08022, 2.98870, 4.05567, 5.27568, 6.63893, 8.12404, 9.65014)

plot f(x) notitle, 0.083567448 + f0(x) / 60000 t "0^{th} state " , 0.334104121 +  f1(x) / 60000 t "1^{th} state " , 0.751095414 + f2(x) / 60000 t "2^{th} state " , 1.33361793 + f3(x) / 60000 t "3^{th} state " , 2.08021998+ f4(x) / 60000 t "4^{th} state " , 2.98869848 + f5(x) / 60000 t "5^{th} state " , 4.05567217 + f6(x) / 60000 t "6^{th} state " , 5.27567530 + f7(x) / 60000 t "7^{th} state " , 6.63893080 + f8(x) / 60000 t "8^{th} state " , 8.12403774 + f9(x) / 60000 t "9^{th} state " , 9.65014458 + f10(x) / 60000 t "10^{th} state "
