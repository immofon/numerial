#include <stdio.h>

double calc_while(int n) {
	double sum = 0;
	int i = 1;
	while(i<=n){
		sum += 1/((double)2*i-1);
		i++;
	}
	return sum;
}

double calc_dowhile(int n) {
	double sum = 0;
	int i = 1;
	do {
		sum += 1/((double)2*i-1);
		i++;
	}while(i <= n);
	return sum;
}

double calc_for(int n) {
	double sum = 0;
	for (int i = 1; i <= n;i++) {
		sum += 1/((double)2*i-1);
	}

	return sum;
}

int main() {
	printf("while: %lf\n",calc_while(100));
	printf("dowhile: %lf\n",calc_dowhile(100));
	printf("for: %lf\n",calc_for(100));

	return 0;
}

