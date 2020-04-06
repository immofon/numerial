# numerial

## How to use it on Linux or Max OS X

1. Copy numerial.c and numerial.h to your project path.
2. Include "numerial.h".
3. Write you own code.
4. Compile and run.

`mat_qr_test.c`
```c
#include <stdio.h>

#include "numerial.h"



#define debug_qr() do {\
	printf(":%d:\n",__LINE__); \
	printf("A=\n"); \
	mat_println(".2",A); \
	printf("Q=\n"); \
	mat_println(".2",Q); \
	printf("R=\n"); \
	mat_println(".2",R); \
	printf("--------------------------------\n");\
} while(0)


int main() {
	double Av[] = {
		2,5,-2,
		1,-1,5,
		4,1,-2,
	};
	mat_t A = new_mat(3,3);
	assert(mat_assign(A,Av));
	mat_t Q = new_mat(3,3);
	mat_t R = new_mat(3,3);
	mat_t b = new_mat_vec(3);
	double bv[] = {1,0,0};
	assert(mat_assign(b,bv));
	mat_t c = new_mat_vec(3);


	assert(mat_reduction_qr(Q,R,A));
	debug_qr();

	assert(mat_transpose(Q,Q));
	mat_println("",Q);

	assert(mat_product(c,Q,b));
	mat_println("",c);

	assert(mat_back_solution(R,b,c));
	mat_println("",b);

	assert(mat_product(c,A,b));
	mat_println("",c);

	assert(mat_inv_qr(R,A));
	mat_println("",R);


	assert(mat_product(Q,A,R));
	mat_println("",Q);


	free_mat(&A);
	free_mat(&Q);
	free_mat(&R);
	return 0;
}

```


```shell
clang mat_qr_test.c numerial.c
./a.out
```

```shell
:36:
A=
[2.00 5.00 -2.00 
 1.00 -1.00 5.00 
 4.00 1.00 -2.00 ]
Q=
[0.44 0.86 0.25 
 0.22 -0.37 0.90 
 0.87 -0.34 -0.35 ]
R=
[4.58 2.84 -1.53 
 0.00 4.35 -2.91 
 0.00 0.00 4.71 ]
--------------------------------
[0.436436 0.218218 0.872872 
 0.864124 -0.371901 -0.339087 
 0.250627 0.902258 -0.350878 ]
[0.436436 
 0.864124 
 0.250627 ]
[-0.031915 
 0.234043 
 0.053191 ]
[1.000000 
 0.000000 
 -0.000000 ]
[-0.031915 0.085106 0.244681 
 0.234043 0.042553 -0.127660 
 0.053191 0.191489 -0.074468 ]
[1.000000 0.000000 0.000000 
 0.000000 1.000000 -0.000000 
 -0.000000 0.000000 1.000000 ]

```
