#include <stdio.h>
#include "numerial.h"

int main() {
	double A_v[] = {
		13,9,9,11,
		9,11,6,9,
		9,6,12,7,
		11,9,7,10,
	};
	mat_t A = new_mat(4,4);
	assert(mat_assign(A,A_v));

	double b_v[] = {42,35,34,37};
	mat_t b = new_mat_vec(4);
	assert(mat_assign(b,b_v));

	mat_t x = new_mat_vec(4);
	iter_conf_t conf;
	conf.tol = 1e-6;
	conf.max_step = 2000;

	assert(mat_solve_iter_steepest_descent(x,A,b,&conf));
	mat_println("",x);

	return 0;
}
