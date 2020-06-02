#include <stdio.h>
#include "numerial.h"

int main() {
	double A_v[] = {
		1.0,0.5,0.4,0.3,0.2,0.1,
		0.5,1.0,0.5,0.4,0.3,0.2,
		0.4,0.5,1.0,0.5,0.4,0.3,
		0.3,0.4,0.5,1.0,0.5,0.4,
		0.2,0.3,0.4,0.5,1.0,0.5,
		0.1,0.2,0.3,0.4,0.5,1.0,
	};
	mat_t A = new_mat(6,6);
	assert(mat_assign(A,A_v));

	printf("A = \n");
	mat_println("4.1",A);

	double xc_v[] = {1,8,0,4,3,5};
	mat_t xc = new_mat_vec(6);
	assert(mat_assign(xc,xc_v));

	printf("xc = \n");
	mat_println("2.0",xc);

	mat_t b = new_mat_vec(6);
	assert(mat_product(b,A,xc));

	printf("b = A * xc = \n");
	mat_println("",b);

	mat_t x = new_mat_vec(6);
	iter_conf_t conf;
	conf.tol = 1e-6;
	conf.max_step = 2000;

	printf("Solve Ax=b\n");
	printf("Method: Steepest Descent:\nx = \n");
	assert(mat_solve_iter_steepest_descent(x,A,b,&conf));
	mat_println("",x);
	printf("used_step: %d\n\n",conf.used_step);

	int i;
	range(i,1,x.m,1) {
		vec_v(x,i) = 0;
	}

	printf("Method: Conjugate Gradient:\nx = \n");
	assert(mat_solve_iter_conjugate_gradient(x,A,b,&conf));
	mat_println("",x);
	printf("used_step: %d\n",conf.used_step);

	pause();
	free_mat(&A);
	free_mat(&b);
	free_mat(&xc);
	free_mat(&x);
	return 0;
}
