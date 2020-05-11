#include <stdio.h>
#include "numerial.h"

int main() {
	mat_t A = new_mat(6,6);
	mat_t L = new_mat(6,6);
	mat_t LT = new_mat(6,6);
	mat_t xc = new_mat_vec(6);
	mat_t x = new_mat_vec(6);
	mat_t y = new_mat_vec(6);
	mat_t b = new_mat_vec(6);

	double A_v[] =  {
		1,1,1,1,1,1,
		1,2,2,2,2,2,
		1,2,3,3,3,3,
		1,2,3,4,4,4,
		1,2,3,4,5,5,
		1,2,3,4,5,6,
	};
	MUST(mat_assign(A,A_v));
	printf("A = \n");
	mat_println(".0",A);

	double xc_v[] = {
		1,8,0,4,3,5,
	};
	MUST(mat_assign(xc,xc_v));
	printf("x_c = \n");
	mat_println(".0",xc);

	MUST(mat_product(b,A,xc));
	printf("b = A * x_c = \n");
	mat_println("",b);

	printf("Cholesky reduction and solve Ax=b got: \n");
	MUST(mat_reduction_llt_cholesky(L,A));
	MUST(mat_transpose(LT,L));

	MUST(mat_back_substitution_L(y,L,b));
	MUST(mat_back_substitution_U(x,LT,y));


	printf("x = \n");
	mat_println("",x);

	pause();

	HANDLE_MUST(ret);

	free_mat(&A);
	free_mat(&L);
	free_mat(&LT);
	free_mat(&xc);
	free_mat(&x);
	free_mat(&y);
	free_mat(&b);

	return ret;
}
