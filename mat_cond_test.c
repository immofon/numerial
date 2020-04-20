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

	printf("cond_1 %lf\n",mat_cond_1(A));
	printf("cond_inf %lf\n",mat_cond_inf(A));
	printf("cond_F %lf\n",mat_cond_F(A));


	printf("cond_1 %lf\n",mat_cond_1(Ainv));
	printf("cond_inf %lf\n",mat_cond_inf(Ainv));
	printf("cond_F %lf\n",mat_cond_F(Ainv));
	return 0;
}
