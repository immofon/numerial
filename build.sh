#!/bin/bash

TOTOL="8"

echo "Build $1" &&
	echo "Target: build/$1" &&
	echo "[1/${TOTOL}]append numerial.h to build/$1 " &&
	cat numerial.h > build/$1 &&
	echo "[2/${TOTOL}]append numerial.c to build/$1 " &&
	cat numerial.c >> build/$1 &&
	echo "[3/${TOTOL}]append $1 to build/$1 " &&
	cat $1 >> build/$1 &&
	echo "[4/${TOTOL}]remove #pragma once in build/$1 " &&
	sed "s/\#pragma once//g" build/$1 > build/$1.tmp &&
	echo "[5/${TOTOL}]remove #include \"numerial.h\" in build/$1 " &&
	sed "s/\#include \"numerial.h\"//g" build/$1.tmp > build/$1 &&
	rm build/$1.tmp &&
	echo "[6/${TOTOL}]compile build/$1 use clang to build/a.out" &&
	clang build/$1 -o build/a.out &&
	echo "[7/${TOTOL}]run build/a.out" &&
	echo "----------------------------------------------------------" &&
	build/a.out &&
	echo "----------------------------------------------------------" &&
	echo "[8/${TOTOL}]remove build/a.out" &&
	rm -f build/a.out
