#!/bin/sh
v="$2"
a="$3"

for i in $(seq 0 100); do
	x0=$(echo "scale = 3;\
		var = $v - (($i * 2 * a(1)) / (0.511974 * $a))^2;
		if (var < 0) var = 0 else var = var; var"\
		| bc -l)
	x1=$(echo "scale = 3;\
		var = $v - ((($i + 1) * 2 * a(1)) / (0.511974 * $a))^2;\
		if (var < 0) var = 0 else var = var; var"\
		| bc -l)
	if [ $(echo "($x0*100)/1" | bc) -eq 0 ] && [ $(echo "($x1*100)/1" | bc) -eq 0 ]; then
		exit
	else
		"$1" "$x0" "$x1" "$i"
	fi

#	"$1" "$v" "$a" "$x0" "$x1" "$i"
done
