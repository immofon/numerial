#include <stdio.h>
#include "numerial.h"

int main() {
	mat_t H = new_mat(4,4);
	mat_t xc = new_mat_vec(4);
	mat_t g = new_mat_vec(4);
	mat_t x = new_mat_vec(4);

	double eps = 1e-6;
	int step = 1000;

	double H_v[] = {
		0.5,0.125,-0.125,0.1,
		0.1,-0.5,0.15,-0.125,
		0.13,0.125,0.5,0.13,
		0.1,0.125,0.13,-0.5,
	};
	assert(mat_assign(H,H_v));

	printf("H = \n");
	mat_println("7.3",H);

	double xc_v[] = {0,4,3,5};
	assert(mat_assign(xc,xc_v));

	printf("x_c = \n");
	mat_println("12.8",xc);

	assert(mat_product(g,H,xc));
	assert(mat_sub(g,xc,g));

	printf("g = x_c - H*x_c = \n");
	mat_println("12.8",g);

	printf("Solve x = H*x + g\n");

	int i,j;
	printf("Simple Iterative Method:\n");
	init_mat(x,i,j,0);
	assert(mat_solve_iter_simple(x,H,g,step,eps));

	printf("x = \n");
	mat_println("12.8",x);

	printf("Seidel Iterative Method:\n");
	init_mat(x,i,j,0);
	assert(mat_solve_iter_simple(x,H,g,step,eps));
	printf("x = \n");
	mat_println("12.8",x);

	pause();
	free_mat(&H);
	free_mat(&xc);
	free_mat(&g);
	free_mat(&x);
	return 0;
}
