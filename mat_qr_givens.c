#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define new(type,size) ((type *) malloc(sizeof(type)*(size)))
#define each(i,n) for (i=0;i<(n);i++)
#define range(i,from,to,delta) for (i=(from);i<=(to);i+=(delta))

typedef struct{
	int m;
	int n;
	double* data;
} mat_t;

// m*n matrix, use it to read and write.
// i,j are both one-based.
#define mat_v(mat,i,j) (mat.data[(((i)-1)*(mat.n)+((j)-1))])

#define mat_each(mat,i,j) range(i,1,(mat).n,1) range(j,1,(mat).n,1)

#define init_mat(mat,i,j,exp) mat_each(mat,i,j) { \
	mat_v(mat,i,j) = (exp); \
}

// MUST call free_mat after using.
mat_t new_mat(int m,int n);
mat_t new_mat_vec(int n); // vertical vector

#define vec_v(mat,i) (mat.data[(i)-1])

void free_mat(mat_t* mat);

// format: "12.8" default if format == NULL
void mat_println(const char* format,mat_t mat);

// R = A*B
// R should not be A or B.
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R,mat_t A,mat_t B);

int mat_assign(mat_t mat,double* data);

// T can be A.
int mat_transpose(mat_t T, mat_t A);
int mat_back_substitution(mat_t A,mat_t x,mat_t b);
int mat_reduction_qr_givens(mat_t Q,mat_t R,mat_t A);
int mat_inv_qr_givens(mat_t T,mat_t A);


void print_double(const char* format,double v) {
	int i;


	size_t len = strlen(format) + 3 + 1; // "%lf" =>  3
	char* fmt = new(char,len);
	each(i,len) {
		fmt[i] = 0;
	}

	strcat(fmt,"%");
	strcat(fmt,format);
	strcat(fmt,"lf");

	printf(fmt,v);

	free(fmt);
}
// MUST call free_mat after using.
mat_t new_mat(int m,int n) {
	assert(m>0);
	assert(n>0);

	mat_t mat;
	mat.m = m;
	mat.n = n;
	mat.data = new(double,m*n);
	return mat;
}
void free_mat(mat_t* mat) {
	if(mat->data == NULL) {
		return;
	}

	free(mat->data);
	mat->data = NULL;
}

// vertical vector
mat_t new_mat_vec(int n)  {
	return new_mat(n,1);
}


int mat_is_same_size(mat_t A,mat_t B) {
	return A.m == B.m && A.n == B.n;
}

int mat_is_same_size_3(mat_t A,mat_t B,mat_t C) {
	return mat_is_same_size(A,B) && mat_is_same_size(B,C);
}

// share same memories.
int mat_is_identical(mat_t A,mat_t B) {
	if(!mat_is_same_size(A,B)) {
		return 0;
	}
	if(A.data == NULL || B.data == NULL) {
		return 0;
	}
	if(A.data == B.data) {
		return 1;
	}
	return 0;
}

// R = A*B
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R,mat_t A,mat_t B) {

	int i,j,k;
	double sum;
	if(mat_is_identical(R,A) || mat_is_identical(R,B)) {
		return 0;
	}
	if (R.m != A.m || R.n != B.n || A.n != B.m) {
		return 0;
	}

	range(i,1,A.m,1) range(j,1,B.n,1) {
		sum = 0;
		range(k,1,A.n,1) {
			sum += mat_v(A,i,k) * mat_v(B,k,j);
		}
		mat_v(R,i,j) = sum;
	}
	return 1;
}

// format: "12.8" default if format == NULL
void mat_println(const char* format,mat_t mat){
	int i,j;

	if (format == NULL ) {
		format = "12.8";
	}

	range(i,1,mat.m,1) {
		if (i == 1) {
			printf("[");
		}else {
			printf(" ");
		}

		range(j,1,mat.n,1) {
			print_double(format,mat_v(mat,i,j));
			printf(" ");
		}

		if (i == mat.m) {
			printf("]");
		}

		printf("\n");
	}
}

int mat_assign(mat_t mat,double* data) {
	assert(mat.data != NULL);
	assert(data != NULL);

	int size = mat.m * mat.n-1;
	int i;
	range(i,0,size,1) {
		mat.data[i] = data[i];
	}
	return 1;
}

int mat_copy(mat_t dst,mat_t src) {
	if (dst.m != dst.m || src.n != src.n) {
		return 0;
	}

	if (dst.data == NULL || src.data == NULL) {
		return 0;
	}

	return mat_assign(dst,src.data);
}

// T can be A.
int mat_transpose(mat_t T, mat_t A) {
	if (T.m != A.n || T.n != A.m) {
		return 0;
	}

	double tmp;
	int i,j;

	if (mat_is_identical(T,A) && A.m == A.n) {
		range(i,1,A.n,1) range(j,i+1,A.n,1) {
			tmp = mat_v(A,i,j);
			mat_v(A,i,j)= mat_v(A,j,i);
			mat_v(A,j,i)= tmp;
		}
		return 1;
	}

	mat_each(T,i,j) {
		mat_v(T,i,j) = mat_v(A,j,i);
	}
	return 1;
}


int mat_back_substitution(mat_t A,mat_t x,mat_t b) {
	if (x.n != 1 || b.n != 1 ||  A.m != A.n || x.m != b.m) {
		return 0;
	}

	int j,k;
	double tmp;

	vec_v(x,A.n) = vec_v(b,A.n) / mat_v(A,A.n,A.n);

	for (k = A.n - 1; k >= 1; k--) {
		tmp = vec_v(b,k);
		range(j,k+1,A.n,1) {
			tmp -= mat_v(A,k,j) * vec_v(x,j);
		}
		vec_v(x,k) = tmp / mat_v(A,k,k);
	}
	return 1;
}

int mat_reduction_qr_givens(mat_t Q,mat_t R,mat_t A) {
	if(R.m != A.m || R.n != A.n) {
		return 0;
	}

	if(R.n != A.m || R.m != R.n) {
		return 0;
	}

	int i,j,k;
	double xi,xj,c,s,u;

	if(!mat_copy(R,A)) {
		return 0;
	}

	init_mat(Q,i,j,i==j?1:0);

	range(k,1,A.n,1) range(i,k+1,A.m,1) {
		xi = mat_v(R,k,k);
		xj = mat_v(R,i,k);
		u = sqrt(xi*xi + xj*xj);

		if (u > 1E-8) {
			c = xi/u;
			s = xj/u;

			range(j,1,A.n,1) {
				xi = mat_v(R,k,j);
				xj = mat_v(R,i,j);
				mat_v(R,k,j) = c*xi + s*xj;
				mat_v(R,i,j) = c*xj - s*xi;

				xi = mat_v(Q,k,j);
				xj = mat_v(Q,i,j);
				mat_v(Q,k,j) = c*xi + s*xj;
				mat_v(Q,i,j) = c*xj - s*xi;
			}
		}
	}

	if(!mat_transpose(Q,Q)) {
		return 0;
	}
	return 1;
}

// T cannot be A.
int mat_inv_qr_givens(mat_t T,mat_t A) {
	if (A.m != A.n || T.m != T.n || A.m != T.m) {
		return 0;
	}

	int i,j,k;
	// use goto to return.
	mat_t Q = new_mat(A.n,A.n);
	mat_t R = new_mat(A.n,A.n);
	mat_t E = new_mat_vec(A.n);
	mat_t b = new_mat_vec(A.n);

	if(!mat_reduction_qr_givens(Q,R,A)) {
		goto ERROR;
	}

	if(!mat_transpose(Q,Q)) {
		goto ERROR;
	}

	range(j,1,A.n,1) {
		range(i,1,A.n,1) {
			vec_v(E,i) = i ==j ? 1:0;
		}

		if(!mat_product(b,Q,E)){
			goto ERROR;
		}

		if(!mat_back_substitution(R,E,b)) {
			goto ERROR;
		}

		range(i,1,A.n,1) {
			mat_v(T,i,j) = vec_v(E,i);
		}
	}
	goto OK;


	int ret;
	if(0==0){
OK:
		ret = 1;
	}else {
ERROR:
		ret = 0;
	}
	free_mat(&Q);
	free_mat(&R);
	free_mat(&E);
	free_mat(&b);
	return ret;
}



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

	free_mat(&A);
	free_mat(&B);
	free_mat(&C);

	return 0;
}
