#!/bin/bash

clang -std=c99 -fPIC  -bundle -undefined dynamic_lookup numerial.c num.c -o num.so
mv num.so testlua/num.so
