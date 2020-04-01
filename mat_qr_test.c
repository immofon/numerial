#include <stdio.h>

#include "numerial.h"



#define debug_qr() do {\
	printf(":%d:\n",__LINE__); \
	printf("A=\n"); \
	mat_println(".2",A); \
	printf("Q=\n"); \
	mat_println(".2",Q); \
	printf("R=\n"); \
	mat_println(".2",R); \
	printf("--------------------------------\n");\
} while(0)


int main() {
	double Av[] = {
		2,5,-2,
		1,-1,5,
		4,1,-2,
	};
	mat_t A = new_mat(3,3);
	assert(mat_assign(A,Av));
	mat_t Q = new_mat(3,3);
	mat_t R = new_mat(3,3);
	mat_t b = new_mat_vec(3);
	double bv[] = {1,0,0};
	assert(mat_assign(b,bv));
	mat_t c = new_mat_vec(3);


	assert(mat_reduction_qr(Q,R,A));
	debug_qr();

	assert(mat_transpose(Q,Q));
	mat_println("",Q);

	assert(mat_product(c,Q,b));
	mat_println("",c);

	assert(mat_back_solution(R,b,c));
	mat_println("",b);

	assert(mat_product(c,A,b));
	mat_println("",c);

	assert(mat_inv_qr(R,A));
	mat_println("",R);


	assert(mat_product(Q,A,R));
	mat_println("",Q);


	free_mat(&A);
	free_mat(&Q);
	free_mat(&R);
	return 0;
}
