#include "laba1.h"

int main()
{
	matrix *A = new matrix;
	A->set_Hilbert(10);
	// A->input("in.txt");
	// A->decompositionLDU();
	// A->forward_subs();
	// A->back_subs();
	A->output("out.txt");
	return 0;
}


