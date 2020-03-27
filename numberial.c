#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "numberial.h"

// Pn(x) = a[0]*x^n + a[1]*x^(n-1) + ... + a[n]
// len(P) = n+1
// n: degree of polynomial
double poly_value(double a[],int n, double x) {
	double v = a[0];

	for(int i = 1; i <= n ; i++) {
		v = x * v + a[i];
	}
	return v;
}


// MUST free after use its return.
double* exp_taylor(int n) {
	// double *a = (double *)malloc(sizeof(double)*(n+1));
	double *a = new(double,n+1);
	a[n] = 1;
	for (int i = n-1; i>=0; i--) {
		a[i] = a[i+1] / (n-i);
	}
	return a;
}


double vec_dot_product(double a[], double b[], int dim) {
	double sum = 0;
	for (int i = 0;i < dim; i++) {
		sum += a[i]*b[i];
	}
	return sum;
}

void print_double(const char* format,double v) {
	int i;


	size_t len = strlen(format) + 3 + 1; // "%lf" =>  3
	char* fmt = new(char,len);
	each(i,len) {
		fmt[i] = 0;
	}

	strcat(fmt,"%");
	strcat(fmt,format);
	strcat(fmt,"lf");

	printf(fmt,v);

	free(fmt);
}
// MUST call free_mat after using.
mat_t new_mat(int m,int n) {
	assert(m>0);
	assert(n>0);

	mat_t mat;
	mat.m = m;
	mat.n = n;
	mat.data = new(double,m*n);
	return mat;
}
void free_mat(mat_t* mat) {
	if(mat->data == NULL) {
		return;
	}

	free(mat->data);
	mat->data = NULL;
}

mat_t new_mat_clone(mat_t src) {
	assert(src.data != NULL);

	mat_t R = new_mat(src.m,src.n);
	assert(mat_copy(R,src));
	return R;
}

int mat_is_same_size(mat_t A,mat_t B) {
	return A.m == B.m && A.n == B.n;
}

int mat_is_same_size_3(mat_t A,mat_t B,mat_t C) {
	return mat_is_same_size(A,B) && mat_is_same_size(B,C);
}

// share same memories.
int mat_is_identical(mat_t A,mat_t B) {
	if(!mat_is_same_size(A,B)) {
		return 0;
	}
	if(A.data == NULL || B.data == NULL) {
		return 0;
	}
	if(A.data == B.data) {
		return 1;
	}
	return 0;
}

// R = A*B
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R,mat_t A,mat_t B) {

	int i,j,k;
	double sum;
	if(mat_is_identical(R,A) || mat_is_identical(R,B)) {
		return 0;
	}
	if (R.m != A.m || R.n != B.n || A.n != B.m) {
		return 0;
	}

	range(i,1,A.m,1) range(j,1,B.n,1) {
		sum = 0;
		range(k,1,A.n,1) {
			sum += mat_v(A,i,k) * mat_v(B,k,j);
		}
		mat_v(R,i,j) = sum;
	}
	return 1;
}


// R = A+B
// R can be A or B.
int mat_add(mat_t R,mat_t A,mat_t B){
	int i,j;
	if(!mat_is_same_size_3(R,A,B)) {
		return 0;
	}

	mat_each(R,i,j) {
		mat_v(R,i,j) = mat_v(A,i,j) + mat_v(B,i,j);
	}
	return 1;
}

// format: "12.8" default if format == NULL
void mat_println(const char* format,mat_t mat){
	int i,j;

	if (format == NULL ) {
		format = "12.8";
	}

	range(i,1,mat.m,1) {
		if (i == 1) {
			printf("[");
		}else {
			printf(" ");
		}

		range(j,1,mat.n,1) {
			print_double(format,mat_v(mat,i,j));
			printf(" ");
		}

		if (i == mat.m) {
			printf("]");
		}

		printf("\n");
	}
}

int mat_assign(mat_t mat,double* data) {
	assert(mat.data != NULL);
	assert(data != NULL);

	int size = mat.m * mat.n-1;
	int i;
	range(i,0,size,1) {
		mat.data[i] = data[i];
	}
	return 1;
}

int mat_copy(mat_t dst,mat_t src) {
	if (dst.m != dst.m || src.n != src.n) {
		return 0;
	}

	if (dst.data == NULL || src.data == NULL) {
		return 0;
	}

	return mat_assign(dst,src.data);
}


