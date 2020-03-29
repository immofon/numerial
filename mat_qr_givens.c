#include <stdio.h>
#include "numberial.h"

int main() {
	double Av[] = {
		2,5,-2,
		1,-1,5,
		4,1,-2,
	};
	mat_t A = new_mat(3,3);
	assert(mat_assign(A,Av));
	mat_t B = new_mat(3,3);
	mat_t C = new_mat(3,3);

	assert(mat_inv_qr_givens(B,A));
	assert(mat_product(C,A,B));

	printf("A=\n");
	mat_println(".2",A);
	printf("B=\n");
	mat_println(".2",B);
	printf("A*B=\n");
	mat_println(".2",C);
	assert(mat_product(C,B,A));
	mat_println(".2",C);

	return 0;
}
