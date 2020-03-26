#pragma once

#include <stdlib.h>
#include <math.h>

#define new(type,size) ((type *) malloc(sizeof(type)*(size)))

#define new_mat(type,m,n) ((type *)malloc(sizeof(type)*(m)*(n)))

// m*n matrix, use it to read and write.
// n: column size
#define mat_v(mat,n,i,j) ((mat)[((i)*(n)+(j))])

// i: int, after each{}, i should equal to n
// n: int, MUST greater than 0.
#define each(i,n) for (i=0;i<(n);i++)


// Pn(x) = a[0]*x^n + a[1]*x^(n-1) + ... + a[n]
// len(P) = n+1
// n: degree of polynomial
double poly_value(double a[],int n, double x);

// MUST free after use its return.
double* exp_taylor(int n);


double vec_dot_product(double a[], double b[], int dim);

// format: "12.8" default if format == 0
void mat_println(const char* format,double* mat,int m,int n);
