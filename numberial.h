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

#define mat_each(mat,i,j) range(i,1,(mat).n,1) range(j,1,(mat).n,1)

#define init_mat(mat,i,j,exp) mat_each(mat,i,j) { \
	mat_v(mat,i,j) = (exp); \
}


// MUST call free_mat after using.
mat_t new_mat(int m,int n);
mat_t new_mat_clone(mat_t src);
mat_t new_mat_vec(int n); // vertical vector

#define vec_v(mat,i) (mat.data[(i)-1])

void free_mat(mat_t* mat);

// format: "12.8" default if format == NULL
void mat_println(const char* format,mat_t mat);

// R = A*B
// R should not be A or B.
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R,mat_t A,mat_t B);

// R = A+B
// R can be A or B.
int mat_add(mat_t R,mat_t A,mat_t B);

int mat_assign(mat_t mat,double* data);
int mat_copy(mat_t dst,mat_t src);
int mat_is_same_size(mat_t A,mat_t B) ;
int mat_is_same_size_3(mat_t A,mat_t B,mat_t C) ;

// share same memories.
int mat_is_identical(mat_t A,mat_t B);

// T can be A.
int mat_transpose(mat_t T, mat_t A);

int mat_reduction_qr(mat_t Q,mat_t R,mat_t A);

int mat_back_solution(mat_t A,mat_t x,mat_t b);

// T cannot be A.
int mat_inv_qr(mat_t T,mat_t A);

int mat_reduction_qr_givens(mat_t Q,mat_t R,mat_t A);
int mat_inv_qr_givens(mat_t T,mat_t A);

// Pn(x) = a[0]*x^n + a[1]*x^(n-1) + ... + a[n]
// len(P) = n+1
// n: degree of polynomial
double poly_value(double a[],int n, double x);

// MUST free after use its return.
double* exp_taylor(int n);


double vec_dot_product(double a[], double b[], int dim);

