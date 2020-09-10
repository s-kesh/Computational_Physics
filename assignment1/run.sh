#!/bin/sh
v="$2"
a="$3"

for i in $(seq 0 20); do
	"$1" $(echo "scale = 3;\
		$v - (($i * 2 * a(1)) / (0.511974 * $a))^2"\
		| bc -l) $(echo "scale = 3;\
		$v - ((($i + 1) * 2 * a(1)) / (0.511974 * $a))^2"\
		| bc -l)  $i
done
