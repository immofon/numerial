#include <stdio.h>
#include "numberial.h"

#define M  2
#define N  2

#define mat_each(mat,i,j) range(i,1,(mat).n,1) range(j,1,(mat).n,1)

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


#define init_mat(mat,i,j,exp) mat_each(mat,i,j) { \
	mat_v(mat,i,j) = (exp); \
}

int mat_assign(mat_t mat,double* data) {
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

	if (A.data == B.data) {
		return 1;
	}

	return 0;
}

// NOTICE: R can be A or B.
int mat_add(mat_t R,mat_t A,mat_t B) {
	int i,j;
	if(!mat_is_same_size_3(R,A,B)) {
		return 0;
	}

	mat_each(R,i,j) {
		mat_v(R,i,j) = mat_v(A,i,j) + mat_v(B,i,j);
	}
	return 1;
}


// NOTICE: R can be A.
int mat_scaler(mat_t R,mat_t A,double a) {
	int i,j;
	if (!mat_is_same_size(R,A)) {
		return 0;
	}

	mat_each(R,i,j) {
		mat_v(R,i,j) = a * mat_v(A,i,j);
	}
	return 1;
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

	mat_t V = new_mat(3,1);
	double V_v[] = {1,2,3};
	assert(mat_assign(V,V_v));
	mat_println("",V);

	free_mat(&A);
	free_mat(&B);
	free_mat(&C);
	free_mat(&R);
	return 0;
}
