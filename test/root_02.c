#include <stdio.h>
#include "numerial.h"

double f(double x) {
	return x*x - x -  2;
}
double f1(double x) {
	return 2*x - 1;
}


int main() {
	double x;
	iter_conf_t conf;
	conf.tol = 1e-6;
	conf.max_step = 1000;

	printf("Method\t\tx\t\tstep\n");

	x = 0.63;
	assert(root_iter_newton(&x,&f,&f1,&conf));
	printf("Newton:\t\t");
	printf("%lf\t%d\n",x,conf.used_step);


	x = 0.63;
	assert(root_iter_newton_down(&x,&f,&f1,0.05,0.01,&conf));
	printf("Newton Down:\t");
	printf("%lf\t%d\n",x,conf.used_step);

	pause();
	return 0;
}

