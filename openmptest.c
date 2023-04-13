#include <stdio.h>

//gcc -fopenmp openmptest.c -o openmptest
int main(void)
{
	#pragma omp parallel
	printf("Hello, world.\n");
	return 0;
}

