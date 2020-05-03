#!/bin/bash

clang -I ./ -O3 numerial.c $1
./a.out
rm -f a.out
