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



// R = A*B
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R,mat_t A,mat_t B) {

	int i,j,k;
	double sum;
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
