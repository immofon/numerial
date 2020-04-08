#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>

// return with resouce free
#define HANDLE_MUST(varname) int varname;\
	if(0==0){ \
		OK: \
		varname = 1; \
	}else { \
		ERROR: \
		varname = 0; \
	}

#define MUST_return_ok() goto OK
#define MUST_return_error() goto ERROR
#define MUST(exp) if(!(exp)) {MUST_return_error();}


#define new(type,size) ((type *) malloc(sizeof(type)*(size)))


// i: int, after each{}, i should equal to n
// n: int, MUST greater than 0.
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

#define mat_each(mat,i,j) range(i,1,(mat).m,1) range(j,1,(mat).n,1)

#define init_mat(mat,i,j,exp) mat_each(mat,i,j) { \
	mat_v(mat,i,j) = (exp); \
}


// MUST call free_mat after using.
// These all will init as 0.
mat_t new_mat(int m,int n);
mat_t new_mat_vec(int n); // vertical vector

mat_t new_mat_clone(mat_t src);

#define vec_v(mat,i) (mat.data[(i)-1])

void free_mat(mat_t* mat);

// format: "12.8" default if format == NULL
void mat_println(const char* format,mat_t mat);

// share same memories.
int mat_is_identical(mat_t A,mat_t B);

// T can be A.
int mat_transpose(mat_t T, mat_t A);

// R = A*B
// R should not be A or B.
// return: 0(error) 1(ok)
// ADVICE: assert(mat_product(R,A,B);
int mat_product(mat_t R,mat_t A,mat_t B);

// R = A+B
// R can be A or B.
int mat_add(mat_t R,mat_t A,mat_t B);

// R = A-B
// R can be A or B.
int mat_sub(mat_t R,mat_t A,mat_t B);

// R = cA
// R can be A.
int mat_scaler(mat_t R,mat_t A,double c);

int mat_assign(mat_t mat,double* data);
int mat_copy(mat_t dst,mat_t src);
int mat_is_same_size(mat_t A,mat_t B) ;
int mat_is_same_size_3(mat_t A,mat_t B,mat_t C) ;

int mat_back_solution(mat_t A,mat_t x,mat_t b);

int mat_reduction_qr(mat_t Q,mat_t R,mat_t A);
int mat_reduction_qr_givens(mat_t Q,mat_t R,mat_t A);
int mat_reduction_qr_household(mat_t Q,mat_t R,mat_t A);

// T cannot be A.
int mat_inv_qr(mat_t T,mat_t A,int (*reduction_qr)(mat_t Q,mat_t R,mat_t A));

int mat_solve_qr(mat_t x,mat_t A,mat_t b, int (*reduction_qr)(mat_t Q,mat_t R,mat_t A));

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

	int i,j;
	mat_t mat;
	mat.m = m;
	mat.n = n;
	mat.data = new(double,m*n);
	init_mat(mat,i,j,0);
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

mat_t new_mat_clone(mat_t src) {
	assert(src.data != NULL);

	mat_t R = new_mat(src.m,src.n);
	assert(mat_copy(R,src));
	return R;
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


// R = A+B
// R can be A or B.
int mat_add(mat_t R,mat_t A,mat_t B){
	int i,j;
	if(!mat_is_same_size_3(R,A,B)) {
		return 0;
	}

	mat_each(R,i,j) {
		mat_v(R,i,j) = mat_v(A,i,j) + mat_v(B,i,j);
	}
	return 1;
}

// R = A+B
// R can be A or B.
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

// R = cA
// R can be A.
int mat_scaler(mat_t R,mat_t A,double c) {
	int i,j;
	if(!mat_is_same_size(R,A)) {
		return 0;
	}

	mat_each(R,i,j) {
		mat_v(R,i,j) = c * mat_v(A,i,j);
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

double _mat_inter_product(mat_t A, mat_t B,int i,int j) {
	assert(A.m == B.m);
	assert(1<=i & i <= A.n);
	assert(1<=j & j <= B.n);

	int k;
	double sum = 0;
	range(k,1,A.m,1) {
		sum += mat_v(A,k,i)*mat_v(B,k,j);
	}
	return sum;
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

int mat_reduction_qr(mat_t Q,mat_t R,mat_t A) {
	if(!mat_is_same_size(Q,A)) {
		return 0;
	}
	if(R.m != A.n || R.n != A.n) {
		return 0;
	}

	int i,j,k,n;
	double tmp;
	mat_t T = new_mat_vec(A.n);

	init_mat(R,i,j,i<=j ? 1:0);

	range(i,1,Q.m,1) {
		mat_v(Q,i,1) = mat_v(A,i,1);
	}

	vec_v(T,1) = _mat_inter_product(Q,Q,1,1);

	range(k,2,Q.n,1) {
		range(i,1,k-1,1) {
			mat_v(R,i,k) = _mat_inter_product(A,Q,k,i) / vec_v(T,i);
		}

		range(i,1,Q.m,1) {
			tmp = mat_v(A,i,k);
			range(j,1,k-1,1) {
				tmp -= mat_v(R,j,k) * mat_v(Q,i,j);
			}
			mat_v(Q,i,k) = tmp;
		}

		vec_v(T,k) = _mat_inter_product(Q,Q,k,k);
	}

	range(i,1,T.m,1) {
		vec_v(T,i) = sqrt(vec_v(T,i));
	}

	mat_each(Q,i,j) {
		mat_v(Q,i,j) /= vec_v(T,j);
	}
	mat_each(R,i,j) {
		mat_v(R,i,j) *= vec_v(T,i);
	}
	return 1;
}


int mat_back_solution(mat_t A,mat_t x,mat_t b) {
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

// T cannot be A.
int mat_inv_qr(mat_t T,mat_t A,int (*reduction_qr)(mat_t Q,mat_t R,mat_t A)){
	if (A.m != A.n || T.m != T.n || A.m != T.m) {
		return 0;
	}

	int i,j,k;
	// use goto to return.
	mat_t Q = new_mat(A.n,A.n);
	mat_t R = new_mat(A.n,A.n);
	mat_t E = new_mat_vec(A.n);
	mat_t b = new_mat_vec(A.n);

	MUST((*reduction_qr)(Q,R,A));
	MUST(mat_transpose(Q,Q));

	range(j,1,A.n,1) {
		range(i,1,A.n,1) {
			vec_v(E,i) = i ==j ? 1:0;
		}

		MUST(mat_product(b,Q,E));
		MUST(mat_back_solution(R,E,b));

		range(i,1,A.n,1) {
			mat_v(T,i,j) = vec_v(E,i);
		}
	}

	HANDLE_MUST(ret);
	free_mat(&Q);
	free_mat(&R);
	free_mat(&E);
	free_mat(&b);
	return ret;
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



int mat_solve_qr(mat_t x,mat_t A,mat_t b, int (*reduction_qr)(mat_t Q,mat_t R,mat_t A)){
	if (reduction_qr == NULL) {
		return 0;
	}

	if (!mat_is_same_size(x,b)) {
		return 0;
	}

	if(A.m != A.n || A.n != b.m) {
		return 0;
	}

	mat_t Q = new_mat(A.m,A.m);
	mat_t R = new_mat(A.m,A.m);
	mat_t b1 = new_mat_vec(b.m);

	MUST((*reduction_qr)(Q,R,A));
	MUST(mat_transpose(Q,Q));
	MUST(mat_product(b1,Q,b));
	MUST(mat_back_solution(R,x,b1));

	HANDLE_MUST(ret);
	free_mat(&Q);
	free_mat(&R);
	free_mat(&b1);
	return ret;
}

int _mat_set_household(mat_t H,int k,mat_t A) {
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
		MUST(_mat_set_household(H,k,R));
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

void pause() {
#ifdef __WINDOWS_
	system("PAUSE");
#endif
#ifdef _WIN32
	system("PAUSE");
#endif
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
		printf("time: %.12lfms\n",(((double)(clock() - t)*1000)/CLOCKS_PER_SEC));
		printf("Ax=b; x=\n");
		mat_println(".12",x);
		printf("\n");

		free_mat(&A);
		free_mat(&b);
		free_mat(&x);
	}

	pause();

	return 0;
}
