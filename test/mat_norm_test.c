#include <stdio.h>
#include "numerial.h"

int main() {
	double A_v[] =  {
		-1,2,
		3,-4,
	};
	mat_t A = new_mat(2,2);
	mat_t Ainv = new_mat(2,2);
	assert(mat_assign(A,A_v));
	assert(mat_inv_qr(Ainv,A,&mat_reduction_qr));

	printf("A=\n");
	mat_println(".0",A);
	printf("inverse of A=\n");
	mat_println("",Ainv);

	printf("norm_1 %lf\n",mat_norm_1(A));
	printf("norm_inf %lf\n",mat_norm_inf(A));
	printf("norm_F %lf\n",mat_norm_F(A));


	printf("norm_1 %lf\n",mat_norm_1(Ainv));
	printf("norm_inf %lf\n",mat_norm_inf(Ainv));
	printf("norm_F %lf\n",mat_norm_F(Ainv));
	return 0;
}
