#include "numerial.h"

int main() {
	double xx[] = {0,0.2,0.4,0.6,0.8,1};
	double yy[] = {0,0.199,0.389,0.565,0.717,0.841};
	mat_t X = new_mat_vec(3);
	mat_t Y = new_mat_vec(3);
	assert(mat_assign(X,xx));
	assert(mat_assign(Y,yy));

	double y;
	double x;
	for (x = -1; x <= 1; x+=0.01) {
		assert(itp_lagrange(&y,X,Y,x));
		printf("%lf\t%lf\n",x,y);
	}
	return 0;
}
