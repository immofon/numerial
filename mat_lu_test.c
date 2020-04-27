#include <stdio.h>
#include "numerial.h"

int main() {
	mat_t L = new_mat(3,3);
	mat_t U = new_mat(3,3);
	mat_t A = new_mat(3,3);
	double A_v[] = {
		2,4,-2,
		1,-1,5,
		4,1,-2,
	};
	mat_t LU = new_mat(3,3);
	assert(mat_assign(A,A_v));
	assert(mat_reduction_lu_crout(L,U,A));
	assert(mat_product(LU,L,U));
	assert(mat_sub(LU,A,LU));

	mat_println("",L);
	mat_println("",U);
	printf("delta norm: %lf\n",mat_norm_inf(LU));

	double det;
	assert(mat_det_lu(&det,A,&mat_reduction_lu_crout));
	printf("det(A)=%lf\n",det);

	return 0;
}
