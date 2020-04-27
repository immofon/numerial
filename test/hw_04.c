#include <stdio.h>
#include "numerial.h"

int main() {
	double A_v[] = {
		2,6,
		2,6.00001,
	};
	mat_t A = new_mat(2,2);
	assert(mat_assign(A,A_v));

	double b_v[] = {8,8.00001};
	mat_t b = new_mat_vec(2);
	assert(mat_assign(b,b_v));

	double c_v[] = {8,7.99998};
	mat_t c = new_mat_vec(2);
	assert(mat_assign(c,c_v));

	// b + db = c <=> db = c - b
	mat_t db = new_mat_vec(2);
	assert(mat_sub(db,c,b));

	mat_t x = new_mat_vec(2);
	mat_t x1 = new_mat_vec(2);
	mat_t dx = new_mat_vec(2);


	assert(mat_solve_qr(x,A,b,&mat_reduction_qr_givens));
	assert(mat_solve_qr(x1,A,c,&mat_reduction_qr_givens));

	// x1 = x + dx <=> dx = x1 - x
	assert(mat_sub(dx,x1,x));

	double norm_x = vec_norm_inf(x);
	double norm_dx = vec_norm_inf(dx);
	double norm_b = vec_norm_inf(b);
	double norm_db = vec_norm_inf(db);

	double condA = 0;
	assert(mat_cond(&condA,A,&mat_norm_inf));

	printf("Solve Ax = b\n");
	printf("A = \n");
	mat_println("",A);
	printf("b = \n");
	mat_println("",b);
	printf("x = \n");
	mat_println("",x);

	printf("\nSolve A(x+dx) = (b+db)\n");
	printf("A = \n");
	mat_println("",A);
	printf("b+db = \n");
	mat_println("",c);
	printf("x+dx = \n");
	mat_println("",x1);

	printf("dx = \n");
	mat_println("",dx);

	printf("\n");

	printf("norm(x) = %lf\n",norm_x);
	printf("norm(dx) = %lf\n",norm_dx);
	printf("norm(b) = %lf\n",norm_b);
	printf("norm(db) = %lf\n",norm_db);

	printf("\n");

	double delta = ((condA * norm_db)/norm_b) - (norm_dx/norm_x);
	printf("((cond(A) * norm(db))/norm(b)) - (norm(dx)/norm(x)) = %lf >= 0\n",delta);
	printf("Which means that\n");
	printf("(norm(dx)/norm(x)) <= (cond(A) * norm(db)/norm(b))\n");

	free_mat(&A);
	free_mat(&b);
	free_mat(&c);
	free_mat(&db);
	free_mat(&x);
	free_mat(&x1);
	free_mat(&dx);

	pause();
	return 0;
}
