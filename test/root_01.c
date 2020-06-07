#include <stdio.h>
#include "numerial.h"

double f1(double x) {
	return x*x - 2;
}

double f2(double x) {
	return sqrt(x+2);
}

double f3(double x) {
	return 1+2/x;
}

double f4(double x) {
	return (x*x+2)/(2*x - 1);
}

void test(const char* name,fn_t fn) {
	printf("Format: %s\n",name);

	double x = 1.5;
	iter_conf_t conf = {0};
	conf.tol = 1e-6;
	conf.max_step = 1000;

	x = 1.5;
	if(root_iter_fixed_point(&x,fn,&conf)) {
		printf("Convergence:\tx = %.6f\n\t\tused_step = %d\n",x,conf.used_step);
	}else {
		printf("Divergence!\n");
	}
	printf("\n");
}

int main() {
	test("x = x*x - 2",f1);
	test("x = sqrt(x + 2)",f2);
	test("x = 1 + 2/x",f3);
	test("x = (x*x + 2) / (2*x - 1)",f4);

	pause();
	return 0;
}
