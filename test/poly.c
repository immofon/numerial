#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerial.h"


int main() {
	double *a;
	int n[] = {1,3,5,10,30,50,90,-1};
	double x[] = {0.1,0.5,1.0,2.0,5.0,-1};

	printf("\t");
	for (int j = 0; x[j] > 0; j++) {
		printf("%lf\t",x[j]);
	}
	printf("\n");

	printf("exp(x)\t");
	for (int j = 0; x[j] > 0; j++) {
		printf("%.10lf\t",exp(x[j]));
	}
	printf("\n");


	for (int i = 0; n[i] > 0; i++) {
		printf("n: %d\t",n[i]);
		for (int j = 0; x[j] > 0; j++) {
			a = exp_taylor(n[i]);
			printf("%.10lf\t",poly_value(a,n[i],x[j]));
			free(a);
		}
		printf("\n");
	}
}
