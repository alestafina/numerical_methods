#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

typedef double real;

struct matrix
{
	real **al = NULL, **au = NULL, *di = NULL;
	real *b = NULL;
	int n = 0, k = 0;

	void testgen();
	void input_testmem();
	void output_inp_data(string name);
	void mult_LDU();

	void input(string name);
	void decompositionLDU();
	void forward_subs();
	void back_subs();
	void output(string name);
};


