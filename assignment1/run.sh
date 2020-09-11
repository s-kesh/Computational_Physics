#!/bin/bash
v="$2"
a="$3"

for i in $(seq 0 100); do
	x0=$(echo "scale = 3;\
		var = (($i * 2 * a(1)) / (0.511974 * $a))^2;
		if (var < 0) var = 0 else var = var; var"\
		| bc -l)
	x1=$(echo "scale = 3;\
		var = ((($i + 1) * 2 * a(1)) / (0.511974 * $a))^2;\
		if (var < 0) var = 0 else var = var; var"\
		| bc -l)

	x0int=$(echo "$x0/1" | bc)
	x1int=$(echo "$x1/1" | bc)
	vint=$(echo "$v/1" | bc)

	if [ $x0int -lt $vint ];	then
		if [ $x1int -gt $vint ]; then
			x1=$v
		fi
		echo "Eigen Value no $i"
		"$1" "$v" "$a" "$x0" "$x1" "$i"
	else
		exit
	fi

done
