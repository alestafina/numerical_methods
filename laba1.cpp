#include "laba1.h"

int main()
{
	matrix *A = new matrix;
	A->input("in.txt");
	A->decompositionLDU();
	A->forward_subs();
	A->back_subs();
	A->output("out.txt");
	return 0;
}


