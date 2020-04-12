#include "numerial.h"
#include <stdio.h>

mat_t new_mat_hilbert(int n) {
	assert(n >= 1);
	mat_t H = new_mat(n,n);
	int i,j;
	range(i,1,n,1) range(j,1,n,1) {
		mat_v(H,i,j) = ((double)1) / ((double)(i+j-1));
	}
	return H;
}

int main() {

	printf("Stage 1: Build Hilbert matrix H\n");
	mat_t H = new_mat_hilbert(6);
	printf("H = \n");
	mat_println("",H);

	printf("\nStage 2: Build vector b\n");
	double Num_v[] = {1,8,0,4,3,5};
	mat_t Num = new_mat_vec(6);
	assert(mat_assign(Num,Num_v));
	printf("Num = \n");
	mat_println("",Num);

	mat_t b = new_mat_vec(6);

	assert(mat_product(b,H,Num));
	printf("b = \n");
	mat_println("",b);


	printf("\nStage 3: Use two methods to solve Hx=b\n");

	mat_t x = new_mat_vec(6);
	printf("Method 1: Givens\n");

	assert(mat_solve_qr(x,H,b,&mat_reduction_qr_givens));
	printf("x = \n");
	mat_println("",x);

	printf("Method 2: Household\n");

	assert(mat_solve_qr(x,H,b,&mat_reduction_qr_household));
	printf("x = \n");
	mat_println("",x);

	return 0;
}

