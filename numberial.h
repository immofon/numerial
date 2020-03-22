#pragma once

#include <stdlib.h>

#define new(type,size) ((type *) malloc(sizeof(type)*(size)))

// Pn(x) = a[0]*x^n + a[1]*x^(n-1) + ... + a[n]
// len(P) = n+1
// n: degree of polynomial
double poly_value(double a[],int n, double x);

// MUST free after use its return.
double* exp_taylor(int n);
