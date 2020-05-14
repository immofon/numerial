#include <stdio.h>
#include "numerial.h"

int main() {
	mat_t A = new_mat(3,3);
	mat_t x = new_mat_vec(3);

	double A_v[] = {
		1,0,1,
		8,0,4,
		3,5,4,
	};
	assert(mat_assign(A,A_v));

	double eigenvalue = 0;
	iter_conf_t conf;
	assert(mat_eigen_power_method(&eigenvalue,x,A,&conf));

	mat_println(".12",x);
	printf("eigen value: %.12lf\n",eigenvalue);
	printf("iter step: %d\n",conf.used_step);

	return 0;
}

