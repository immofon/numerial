#include <stdio.h>
#include "numerial.h"

int mat_solve_iter_sor_m(mat_t x,mat_t A,mat_t b,double factor,iter_conf_t *conf) {
	mat_t D = new_mat(A.m,A.n);
	mat_t L = new_mat(A.m,A.n);
	mat_t U = new_mat(A.m,A.n);
	mat_t B = new_mat(A.m,A.n);
	mat_t M = new_mat(A.m,A.n);
	mat_t Mb = new_mat_vec(b.m);
	mat_t BM = new_mat(A.m,A.n);

	MUST(A.m == A.n);
	MUST(A.n == b.m);
	MUST(b.n == 1);
	MUST(mat_is_same_size(x,b));

	int i,j;
	int n = A.m;
	// init
	range(i,1,n,1) {
		mat_v(D,i,i) = mat_v(A,i,i);
	}
	range(i,1,n,1) range(j,i+1,n,1) {
		mat_v(U,i,j) = -mat_v(A,i,j);
	}
	range(j,1,n,1) range(i,j+1,n,1) {
		mat_v(L,i,j) = -factor * mat_v(A,i,j);
	}

	MUST(mat_sub(BM,D,L));
	MUST(mat_inv_qr(B,BM,&mat_reduction_qr_household));

	MUST(mat_scaler(D,D,1-factor));
	MUST(mat_scaler(U,U,factor));
	MUST(mat_add(M,D,U));
	MUST(mat_product(BM,B,M));

	MUST(mat_scaler(B,B,factor));
	MUST(mat_product(Mb,B,b));

	MUST(mat_solve_iter_simple(x,BM,Mb,conf));

	HANDLE_MUST(ret);
	free_mat(&D);
	free_mat(&L);
	free_mat(&U);
	free_mat(&B);
	free_mat(&M);
	free_mat(&BM);
	free_mat(&Mb);
	return ret;
}

int main() {
	double A_v[] = {
		8,-1,1,
		2,10,-1,
		1,1,-5,
	};
	mat_t A = new_mat(3,3);
	assert(mat_assign(A,A_v));

	printf("A = \n");
	mat_println("4.0",A);

	double b_v[] = {1,4,3};
	mat_t b = new_mat_vec(3);
	assert(mat_assign(b,b_v));

	printf("b = \n");
	mat_println("2.0",b);

	double x_v[] = {0.125,0.4,-0.6};
	mat_t x = new_mat_vec(3);
	assert(mat_assign(x,x_v));

	printf("x0 = \n");
	mat_println("5.3",x);

	iter_conf_t conf;
	conf.max_step = 1000;
	conf.tol = 1e-6;
	printf("Solve Ax=b\n");
	printf("SOR:\n");
	assert(mat_solve_iter_sor(x,A,b,1.005,&conf));
	printf("x = \n");
	mat_println("10.7",x);
	printf("used step: %d\n\n",conf.used_step);

	printf("Solve Ax=b\n");
	printf("SOR_M:\n");
	assert(mat_assign(x,x_v));
	assert(mat_solve_iter_sor_m(x,A,b,1.005,&conf));
	printf("x = \n");
	mat_println("10.7",x);
	printf("used step: %d\n",conf.used_step);

	pause();
	free_mat(&A);
	free_mat(&b);
	free_mat(&x);
	return 0;
}
