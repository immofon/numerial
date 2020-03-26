#include <string.h>
#include <stdio.h>
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


void mat_println(const char* format,double* mat,int m,int n) {
	int i,j;

	if (format == 0 ) {
		format = "12.8";
	}

	each(i,m) {
		if (i == 0) {
			printf("[");
		}else {
			printf(" ");
		}

		each(j,n) {
			print_double(format,mat_v(mat,n,i,j));
			printf(" ");
		}

		if (i == (m-1)) {
			printf("]");
		}

		printf("\n");
	}
}
