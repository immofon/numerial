#!/bin/bash

clang numberial.c $1
./a.out
rm -f a.out
