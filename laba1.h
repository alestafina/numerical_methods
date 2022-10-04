#ifndef LABA1_H_
#define LABA1_H_

#include <iostream>
#include <fstream>
#include <string.h>

#define DDF

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

struct matrix
{
	// al - массив для хранения нижнего треугольника матрицы,
	// au - массив для хранения верхнего треугольника матрицы,
	// di - массив для хранения диагональных элементов,
	// b - вектор правых частей, n - размерность матрицы А,
	// k - ширина ленты

	real **al = NULL, **au = NULL, *di = NULL;
	real *b = NULL;
	int n = 0, k = 0;

	void input(string name);
	void decompositionLDU();
	void forward_subs();
	void back_subs();
	void output(string name);
};

#endif // LABA1_H_ 
