#include <stdio.h>
#include <stdlib.h>

int TestAppLAPACK(int argc, char *argv[]);
int TestAppCCS   (int argc, char *argv[]);
int TestAppSLEPC (int argc, char *argv[]);

int main(int argc, char *argv[]) 
{
	//TestAppLAPACK(argc, argv);
	TestAppCCS(argc, argv);
	//TestAppSLEPC(argc, argv);
	return 0;
}

