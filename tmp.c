#include <stdio.h>

double dot_product(double a[], double b[], int dim) {
	double sum = 0;
	for (int i = 0;i < dim; i++) {
		sum += a[i]*b[i];
	}
	return sum;
}

int main() {
	double a[] = {1,2,3};
	double b[] = {2,3,4};
	printf("((1,2,3),(2,3,4))=%lf\n",dot_product(a,b,3));
	// OUTPUT
	// ((1,2,3),(2,3,4))=20.000000
	return 0;
}
