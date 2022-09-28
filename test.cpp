#include "laba1.h"

void matrix::testgen() {
	srand(time_t(NULL));

   for (int i = 0; i < k; i++)
      for (int j = 0; j < k - i; j++) {
         al[i][j] = 0;
         au[i][j] = 0;
      }
   
   for (int i = 1; i < n; i++)
      for (int h = min(i, k), j = k - h; j < k; j++) {
         al[i][j] = rand() % 150;
         au[i][j] = rand() % 150;
      }
   
   for (int i = 0; i < n; i++)
      di[i] = 1.0 + rand() % 10;
   
   for (int i = 0; i < n; i++)
      b[i] = rand() % 50;
         
}

void matrix::input_testmem()
{
	srand(time_t(NULL));
	while (n % 2 != 1 || k % 2 != 1 || k > n) {
		n = 3 + rand() % 10;
		k = 3 + rand() % (n - 3 + 1);
	}

	k = (k - 1) / 2;

	al = new real * [n];
	au = new real * [n];
	di = new real[n];
	b = new real[n];

	for (int i = 0; i < n; i++)
	{
		al[i] = new real[k];
		au[i] = new real[k];
	}
}

void matrix::output_inp_data(string name) {
	ofstream fout(name);

	fout << n << " " << 2 * k + 1 << endl << endl;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < k; j++)
			fout << al[i][j] << " ";
		fout << endl;
	}
	fout << endl;

	for (int i = 0; i < n; i++)
		fout << di[i] << " ";
	fout << endl << endl;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < k; j++)
			fout << au[i][j] << " ";
		fout << endl;
	}
	fout << endl;

	for (int i = 0; i < n; i++)
		fout << b[i] << " ";
}

void matrix::mult_LDU() {

	// Перемножение нижнетреугольной матрицы 
	// и диагональной матрицы (L * D)

	for (int i = k - 1, t = 0; i >= 0; i--, t++)
		for (int j = 1 + t, h = 0; j < n; j++, h++)
			al[j][i] *= di[h];

	// Перемножение новой нижнетреугольной
	// и верхнетреугольной матриц ((LD) * U)

	real sum = 0.0;

	for (int i = 0; i < n; i++)
		for (int j = 0; j < k; j++) {
			
		}
}
