#include "laba1.h"

int main()
{
	matrix *A = new matrix;
	dense_matrix *D = new dense_matrix;
	// A->set_Hilbert(13);
	D->input_DM("in.txt");
	D->Gauss_method();
	// D->band_to_dense(A);
	// A->mult_mx_v();
	// A->decompositionLDU();
	// A->forward_subs();
	// A->back_subs();
	D->output_DM("out.txt");
	return 0;
}


