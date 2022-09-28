#include "laba1.h"

int main()
{
	matrix *A = new matrix;
	A->input_testmem();
	A->testgen();
	A->output("expect.txt");
	A->mult_LDU();
	A->output_inp_data("in.txt");
	return 0;
}


