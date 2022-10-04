#include "laba1.h"

void matrix::input(string name)
{
	ifstream fin(name);
	fin >> n >> k;

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

	for (int i = 0; i < n; i++)
		for (int j = 0; j < k; j++)
			fin >> al[i][j];

	for (int i = 0; i < n; i++)
		fin >> di[i];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < k; j++)
			fin >> au[i][j];

	for (int i = 0; i < n; i++)
		fin >> b[i];

}

void matrix::decompositionLDU()
{
	real sumL, sumU, sumD;

	for (int i = 1; i < n; i++)
	{
		sumD = 0.0;
		for (int j = 0; j < k; j++)
		{
			sumL = 0.0;
			sumU = 0.0;
			for (int h = 1, t = i - k, p = k - j; h <= j && t <= i - j + 1 && p < k; h++, t++, p++)
			{
				sumL += al[i][h - 1] * di[t] * au[abs(i - k + j)][p];
				sumU += au[i][h - 1] * di[t] * al[abs(i - k + j)][p];
			}
			al[i][j] = (al[i][j] - sumL) / di[max(0, i + j - k)];
			au[i][j] = (au[i][j] - sumU) / di[max(0, i + j - k)];
		}
		for (int j = 0, h = k; j < k && h > 0; j++, h--)
		{
			sumD += al[i][j] * di[abs(i - h)] * au[i][j];
		}
		di[i] -= sumD;
	}
}

void matrix::forward_subs()
{
	// Реализация прямого хода: решение Ly = b, 
	// где L - нижнетреугольная матрица с 1 на диагонали,
	// y = DUx, записываем результат в b.
	// Затем решаем Dz = y, D - диагональная матрица, z = Ux, 
	// результат записываем в y (то есть в вектор b).

	real_sum sum;

	for (int i = 1; i < n; i++)
	{
		sum = 0.0;
		for (int j = 0, h = k; j < k && h > 0; j++, h--)
			sum += al[i][j] * b[i - h];
		b[i] -= sum;
	}
	for (int i = 0; i < n; i++)
		b[i] = b[i] / di[i];
}

void matrix::back_subs()
{
	// Реализация обратного хода: решаем систему Ux = z,
	// U - верхнетреугольная матрица с 1 на диагонали,
	// x - решение СЛАУ, ответ записываем в z (то есть опять в вектор b).

	real_sum sum;

	for (int i = n - 2; i >= 0; i--)
	{
		sum = 0.0;
		for (int j = k - 1, h = i + 1; j >= 0 && h < n; j--, h++)
			sum += au[h][j] * b[h];
		b[i] -= sum;
	}
}

void matrix::output(string name)
{
	ofstream fout(name);

	// fout << scientific;
	string s = typeid(b[0]).name();
	fout << fixed;
	if (s == "float")
		fout.precision(7);
	else
		fout.precision(15);

	fout << "al = { " << endl;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < k; j++) {
			fout.width(15);
			fout << al[i][j] << " ";
		}
		if (i != n - 1)
			fout << endl;
	}
	fout << "}" << endl << endl << "di = { ";

	for (int i = 0; i < n; i++)
		fout << di[i] << " ";
	fout << "}" << endl << endl << "au = { " << endl;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < k; j++) {
			fout.width(15);
			fout << au[i][j] << " ";
		}
		if (i != n - 1)
			fout << endl;
	}
	fout << "}" << endl << endl;

	fout << "x = { " << endl;
	for (int i = 0; i < n; i++)
		fout << b[i] << endl;
	fout << "}" << endl;
}

void matrix::set_Hilbert(int size) {
	n = size;
	k = n - 1;

	real *x = new real[n];

	al = new real * [n];
	au = new real * [n];
	di = new real[n];
	b = new real[n];

	for (int i = 0; i < n; i++)
	{
		al[i] = new real[k];
		au[i] = new real[k];
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0, h = 1; j < k; j++) {
			al[i][j] = (i + j < k) ? 0 : 1.0 / (i + (real)h);
			au[i][j] = (i + j < k) ? 0 : 1.0 / (i + (real)h);
			h += al[i][j] != 0 ? 1 : 0;
		}
		di[i] = 1.0 / (2.0 * i + 1);
	}

	for (int i = 0; i < n; i++)
		x[i] = (real)(i + 1);
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < k; j++) {

		}
	}
}