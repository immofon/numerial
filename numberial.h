#pragma once

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define new(type,size) ((type *) malloc(sizeof(type)*(size)))


// i: int, after each{}, i should equal to n
// n: int, MUST greater than 0.
#define each(i,n) for (i=0;i<(n);i++)

#define range(i,from,to,delta) for (i=(from);i<=(to);i+=(delta))


typedef struct{
	int m;
	int n;
	double* data;
} mat_t;

// m*n matrix, use it to read and write.
// i,j are both one-based.
#define mat_v(mat,i,j) (mat.data[(((i)-1)*(mat.n)+((j)-1))])

// MUST call free_mat after using.
mat_t new_mat(int m,int n);
void free_mat(mat_t* mat);

// format: "12.8" default if format == NULL
void mat_println(const char* format,mat_t mat);

// R = A*B
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R,mat_t A,mat_t B);

// Pn(x) = a[0]*x^n + a[1]*x^(n-1) + ... + a[n]
// len(P) = n+1
// n: degree of polynomial
double poly_value(double a[],int n, double x);

// MUST free after use its return.
double* exp_taylor(int n);


double vec_dot_product(double a[], double b[], int dim);
