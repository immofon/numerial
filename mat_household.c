#include <stdio.h>
#include "numerial.h"

mat_t new_mat_tri(int n) {
	assert(n > 0);

	mat_t T = new_mat(n,n);
	int i;
	range(i,1,n,1) {
		mat_v(T,i,i) = -4;
	}

	range(i,1,n-1,1) {
		mat_v(T,i,i+1) = 1;
		mat_v(T,i+1,i) = 1;
	}

	return T;
}

int main() {

	int n = 20;
	clock_t t;

	mat_t A = new_mat_tri(n);
	mat_println(".0",A);
	mat_t b = new_mat_vec(n);
	mat_t x = new_mat_vec(n);
	int i;

	vec_v(b,1) = -27;
	range(i,2,n,1) {
		vec_v(b,i) = -15;
	}


	t = clock();
	assert(mat_solve_qr(x,A,b,&mat_reduction_qr));
	mat_println(".12",x);

	assert(mat_product(b,A,x));
	mat_println("",b);

	printf("time: %fms\n",(((double)(clock() - t)*1000)/CLOCKS_PER_SEC));

	return 0;
}
