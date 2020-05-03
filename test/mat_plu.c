#include "numerial.h"

int main() {
	double A_v[] = {
		2,4,-2,
		1,-1,5,
		4,1,-2,
	};
	mat_t A = new_mat(3,3);
	assert(mat_assign(A,A_v));
	mat_t P = new_mat(3,3);
	mat_t L = new_mat(3,3);
	mat_t U = new_mat(3,3);
	assert(mat_reduction_plu_doolittle(P,L,U,A));
	mat_println(".0",A);
	mat_println("",P);
	mat_println("",L);
	mat_println("",U);
	return 0;
}
