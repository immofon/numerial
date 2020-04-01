#!/bin/bash

clang numerial.c $1
./a.out
rm -f a.out
