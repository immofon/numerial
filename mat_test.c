#include <stdio.h>
#include "numberial.h"

#define M  2
#define N  2

double mat_norm_1(double* mat,int m,int n) {
	double norm = 0.0;
	int i,j;

	each(i,m) each(j,n) {
		norm += fabs(mat_v(mat,n,i,j));
	}
	return norm;
}

double mat_norm_2(double* mat,int m,int n) {
	double norm = 0.0;
	int i,j;

	each(i,m) each(j,n) {
		norm += mat_v(mat,n,i,j)*mat_v(mat,n,i,j);
	}
	return sqrt(norm);
}

// return: MUST free after using.
double* mat_product(double* A,double* B,int m,int n, int p) {
	double* R = new_mat(double,m,n);
	int i,j,k;
	double sum;

	each(i,m) each(j,p) {
		sum = 0;
		each(k,n) {
			sum += mat_v(A,n,i,k) * mat_v(B,p,k,j);
		}
		mat_v(R,p,i,j) = sum;
	}

	return R;
}

#define init_mat(mat,m,n,i,j,exp) each(i,(m)) each(j,(n)) { \
	mat_v((mat),(n),(i),(j)) = (exp); \
}

int main() {
	int m = M;
	int n = N;
	double* A = new_mat(double,m,n);
	int i,j;
	each(i,m) each(j,n) {
		mat_v(A,n,i,j) = i+j;
	}

	mat_println("",A,m,n);

	printf("%lf\n",mat_norm_2(A,m,n));

	free(A);

	double* B = new_mat(double,2,2);
	double* C = new_mat(double,2,2);

	init_mat(B,2,2,i,j,i+j);
	init_mat(C,2,2,i,j,i+j);

	mat_println("",B,2,2);
	mat_println("",C,2,2);

	double* R = mat_product(B,C,2,2,2);
	mat_println("",R,2,2);

	free(B);
	free(C);
	free(R);
	return 0;
}
