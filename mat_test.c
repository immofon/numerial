#include <stdio.h>
#include "numberial.h"

#define M  2
#define N  2


double mat_norm_1(mat_t mat) {
	double norm = 0.0;
	int i,j;

	mat_each(mat,i,j) {
		norm += fabs(mat_v(mat,i,j));
	}
	return norm;
}

double mat_norm_2(mat_t mat) {
	double norm = 0.0;
	int i,j;

	mat_each(mat,i,j) {
		norm += mat_v(mat,i,j)*mat_v(mat,i,j);
	}
	return sqrt(norm);
}


int main() {
	int m = M;
	int n = N;
	mat_t A = new_mat(m,n);
	mat_t B = new_mat(2,2);
	mat_t C = new_mat(2,2);
	mat_t R = new_mat(2,2);
	int i,j;

	mat_each(A,i,j) {
		mat_v(A,i,j) = i+j;
	}
	mat_println("",A);

	init_mat(B,i,j,10);
	mat_v(B,1,2) = 1;
	mat_v(B,2,1) = 1;

	mat_println("",B);

	assert(mat_product(R,A,B));
	mat_println("",R);

	double C_v[] = {
		1,2,
		3,3.14,
	};
	assert(mat_assign(C,C_v));
	mat_println("",C);

	mat_t E = new_mat_clone(C);
	mat_println("",E);
	free_mat(&E);

	assert(mat_add(R,A,B));
	mat_println("",R);

	mat_t V = new_mat(2,1);
	double V_v[] = {1,3};
	assert(mat_assign(V,V_v));
	mat_println("",V);

	mat_t W = new_mat(2,1);
	assert(mat_product(W,B,V));
	mat_println("",W);

	mat_t I = new_mat(10,10);
	init_mat(I,i,j,i==j?1:0);
	mat_println(".0",I);

	free_mat(&A);
	free_mat(&B);
	free_mat(&C);
	free_mat(&R);
	free_mat(&V);
	free_mat(&W);
	return 0;
}
