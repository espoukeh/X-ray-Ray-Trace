#include <stdio.h>
#include "vector.h"

int main(int argc, char *argv[])
{
	char test1[100];
	int i = 4;
	double f = 1.2345678901234567890123;
	sprintf(test1,"this is a test %%1.%ile\n",i);
	printf("%s",test1);
	printf(test1,f);
	return 0;
} 
