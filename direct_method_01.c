#include <stdio.h>
#include "numerial.h"

mat_t new_mat_tri(int n) {
	assert(n>1);
	mat_t A = new_mat(n,n);
	int i,j;
	mat_each(A,i,j) {
		if (i == j) {
			mat_v(A,i,i) = -4;
		}
		if (i+1 == j) {
			mat_v(A,i,j) = 1;
		}
		if (i == j+1) {
			mat_v(A,i,j) = 1;
		}
	}
	return A;
}

int main() {
	double det;
	mat_t A = new_mat_tri(5);
	printf("A = \n");
	mat_println(".0",A);

	printf("Doolittle:\n");
	assert(mat_det_lu(&det,A,&mat_reduction_lu_doolittle));
	printf("det(A)=%lf\n",det);

	printf("Crout:\n");
	assert(mat_det_lu(&det,A,&mat_reduction_lu_crout));
	printf("det(A)=%lf\n",det);

	pause();

	free_mat(&A);
	return 0;
}
