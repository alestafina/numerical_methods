#ifndef LABA1_H_
#define LABA1_H_

#include <iostream>
#include <fstream>
#include <string.h>

#define DD
#ifdef DDF
typedef double real_sum;
typedef float real;
#endif // DDF

#ifdef DF
typedef float real, real_sum;
#endif // DF

#ifdef DD 
typedef double real, real_sum;
#endif // DD

using namespace std;


// com

struct matrix
{
	// al - ������ ��� �������� ������� ������������ �������,
	// au - ������ ��� �������� �������� ������������ �������,
	// di - ������ ��� �������� ������������ ���������,
	// b - ������ ������ ������, n - ����������� ������� �,
	// k - ������ �����, x - ������-������� ���� (����� ��� ������������)

	real **al = NULL, **au = NULL, *di = NULL;
	real *b = NULL, *x = NULL;
	int n = 0, k = 0;

	void input(string name);
	void decompositionLDU();
	void forward_subs();
	void back_subs();
	void output(string name);
	void mult_mx_v();
	void set_Hilbert(int k);
};

struct dense_matrix 
{
	real **a = NULL;
	real *b = NULL, *x = NULL;
	int n;

	void band_to_dense(matrix *m);
	void input_DM(string name);
	void Gauss_method();
	void output_DM(string name);
};

#endif // LABA1_H_ 
