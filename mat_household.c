#include <stdio.h>
#include "numerial.h"

mat_t new_mat_tri(int n) {
	assert(n > 0);

	mat_t T = new_mat(n,n);
	int i,j;
	range(i,1,n,1) {
		mat_v(T,i,i) = -4;
	}

	range(i,1,n-1,1) {
		mat_v(T,i,i+1) = 1;
		mat_v(T,i+1,i) = 1;
	}

	return T;
}





int main() {
	int n;
	int N[] = {5,10,20,-1};
	clock_t t;
	int k = 0;
	int i;

	mat_t A,b,x;
	for (k = 0; N[k] > 0;k++) {
		n = N[k];
		A = new_mat_tri(n);
		b = new_mat_vec(n);
		x = new_mat_vec(n);
		vec_v(b,1) = -27;
		range(i,2,n,1) {
			vec_v(b,i) = -15;
		}

		printf("n=%d\n",n);

		t = clock();
		assert(mat_solve_qr(x,A,b,&mat_reduction_qr));
		printf("Ax=b; x=\n");
		mat_println(".12",x);

		printf("time: %fms\n\n",(((double)(clock() - t)*1000)/CLOCKS_PER_SEC));
		free_mat(&A);
		free_mat(&b);
		free_mat(&x);
	}


	return 0;
}
