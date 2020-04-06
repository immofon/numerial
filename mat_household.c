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

int mat_sub(mat_t R,mat_t A,mat_t B) {
	if(!mat_is_same_size_3(R,A,B)) {
		return 0;
	}

	int i,j;
	mat_each(R,i,j) {
		mat_v(R,i,j) = (mat_v(A,i,j)) - (mat_v(B,i,j));
	}
	return 1;
}

int mat_set_household(mat_t H,int k,mat_t A) {
	if(H.m != H.n) {
		return 0;
	}
	if(H.n != A.m) {
		return 0;
	}
	if(k<1 || k> A.n) {
		return 0;
	}


	mat_t x = new_mat_vec(A.m);
	mat_t y = new_mat_vec(A.m);
	mat_t u = new_mat_vec(A.m);
	mat_t uT = new_mat(1,A.m);
	int i,j=0;
	double alpha=0;
	double sum=0;
	double maxx=0;
	double t=0;

	// set x
	range(i,1,A.m,1) {
		vec_v(x,i) = mat_v(A,i,k);
	}

	// set maxx;
	maxx = fabs(vec_v(x,k));
	range(i,k+1,x.m,1) {
		t = fabs(vec_v(x,i));
		if (maxx < t){
			maxx = t;
		}
	}

	// set alpha
	//TODO lambda
	alpha = 0;
	range(j,k,A.m,1) {
		t = vec_v(x,j)/maxx;
		t *= t;
		alpha += t;
	}
	alpha = maxx*sqrt(alpha);
	if(vec_v(x,k) < -1e-15) {
		alpha = -alpha;
	}

	// set y
	range(i,1,k-1,1) {
		vec_v(y,i) = vec_v(x,i);
	}
	vec_v(y,k) = alpha;
	range(i,k+1,y.m,1) {
		vec_v(y,i) = 0;
	}

	// set sum
	sum = -alpha*(alpha - vec_v(x,k));

	// set u
	MUST(mat_sub(u,x,y));

	MUST(mat_transpose(uT,u));
	if(fabs(sum) < 1e-15) {
		init_mat(H,i,j,i==j?1:0);
		MUST_return_ok();
	}

	MUST(mat_product(H,u,uT));

	sum = 1/sum;

	MUST(mat_scaler(H,H,sum));
	range(i,1,H.m,1) {
		mat_v(H,i,i) += 1;
	}

	HANDLE_MUST(ret);
	free_mat(&x);
	free_mat(&y);
	free_mat(&u);
	free_mat(&uT);
	return ret;
}

int mat_reduction_qr_household(mat_t Q,mat_t R,mat_t A) {
	if(!mat_is_same_size(R,A)) {
		return 0;
	}
	if(Q.m != A.m || Q.n != A.m) {
		return 0;
	}

	int i,j,k;
	mat_t H = new_mat(A.m,A.m);
	mat_t Tq = new_mat(A.m,A.m);
	mat_t Tr = new_mat(A.m,A.n);

	init_mat(Q,i,j,i==j?1:0);
	init_mat(R,i,j,mat_v(A,i,j));

	range(k,1,A.n-1,1) {
		MUST(mat_set_household(H,k,R));
		MUST(mat_product(Tq,H,Q));
		MUST(mat_copy(Q,Tq));
		MUST(mat_product(Tr,H,R));
		MUST(mat_copy(R,Tr));
	}

	MUST(mat_transpose(Q,Q));

	HANDLE_MUST(ret);
	free_mat(&H);
	free_mat(&Tq);
	free_mat(&Tr);
	return ret;
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
