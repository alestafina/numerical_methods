#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

typedef double real;

struct matrix
{
	real **al = NULL, **au = NULL, *di = NULL;
	real *b = NULL;
	int n = 0, k = 0;

	void input(string name);
	void decompositionLDU();
	void forward_subs();
	void back_subs();
	void output(string name);
};


