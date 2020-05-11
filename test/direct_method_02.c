#include <stdio.h>
#include "numerial.h"


int main() {
	mat_t A = new_mat(3,3);
	double A_v[] = {
		2,4,-2,
		1,-1,5,
		4,1,-2,
	};
	assert(mat_assign(A,A_v));

	mat_t P = new_mat(3,3);
	mat_t L = new_mat(3,3);
	mat_t U = new_mat(3,3);

	assert(mat_reduction_plu_doolittle(P,L,U,A));

	int k,i,j;
	mat_t b = new_mat_vec(3);
	mat_t Pb = new_mat_vec(3);
	mat_t x = new_mat_vec(3);
	mat_t y = new_mat_vec(3);
	mat_t invA = new_mat(3,3);
	range(k,1,A.m,1) {
		init_mat(b,i,j,i==k?1:0);
		assert(mat_product(Pb,P,b));
		assert(mat_back_substitution_L(y,L,Pb));
		assert(mat_back_substitution_U(x,U,y));
		range(j,1,x.m,1) {
			mat_v(invA,j,k) = vec_v(x,j);
		}
	}

	printf("A =\n");
	mat_println(".0",A);

	printf("inv(A) =\n");
	mat_println("",invA);

	mat_t I = new_mat(3,3);
	assert(mat_product(I,A,invA));
	printf("A * inv(A) =\n");
	mat_println("",I);

	pause();
	free_mat(&A);
	free_mat(&P);
	free_mat(&L);
	free_mat(&U);
	free_mat(&b);
	free_mat(&Pb);
	free_mat(&x);
	free_mat(&y);
	free_mat(&invA);
	free_mat(&I);
	return 0;
}
