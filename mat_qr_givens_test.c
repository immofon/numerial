#include <stdio.h>
#include "numerial.h"

int main() {
	mat_t A = new_mat(4,3);
	mat_t Q = new_mat(4,4);
	mat_t R = new_mat(4,3);

	double A_v[] = {
		2,4,-2,
		1,-1,5,
		4,1,-2,
		3,5,-7,
	};
	assert(mat_assign(A,A_v));

	assert(mat_reduction_qr_givens(Q,R,A));

	printf("A =\n");
	mat_println(".0",A);

	printf("Q =\n");
	mat_println("",Q);

	printf("R =\n");
	mat_println("",R);

	assert(mat_product(A,Q,R));
	printf("A =\n");
	mat_println("",A);

	mat_t I = new_mat(4,4);
	mat_t QT = new_mat(4,4);

	printf("Make sure Q is orch\n");
	assert(mat_transpose(QT,Q));
	assert(mat_product(I,Q,QT));
	printf("I =\n");
	mat_println("",I);

	return 0;
}
