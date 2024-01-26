

#include "Pippo.h"

#include "Matrix.h"

#include <iostream>


//import Matrix;

//#include <math.h>
//#include <complex.h>
//#include <tgmath.h>

#pragma region USING

// Namespace necessari (evitare: using namespace std;)
using std::cout;
using std::endl;
using std::getchar;
using std::operator<<;
using std::ostream;
using std::complex;
using namespace matrix;
using namespace matrix::linearsys;

using std::string;
#pragma endregion

// Prototyping
void stop_msg();
static int ramdomint(int i, int j);
static int ramdomint10(int i, int j);
static std::complex<double> cpxrc(int r, int c);
static float ramdomfloat10(int i, int j);
static double ramdomdouble20(int i, int j);

int main()
	{

	srand((unsigned)time(NULL));


	cout << "\n\nMatrix test\n" << endl;

	try
		{
		Matrix<Pippo> P1;
		cout << "P1:\n" << P1.to_string() << endl;
		cout << "P1:\n" << P1 << endl;
		Matrix<Pippo> *P2;
		P2 = new Matrix<Pippo>(2, 3);
		(*P2)(0, 1) = Pippo(1);
		(*P2)(1, 2) = Pippo(-5);

		Pippo xpip = (*P2)(1, 2);
		cout << "xpip:\t" << xpip << endl;

		(*P2)(1, 0) = Pippo(-8);
		cout << "P2:\n" << P2->to_string() << "\tRighe: " << P2->rows() << "\tColonne: " << P2->cols() << endl;
		cout << "P2:\n" << P2->to_string(MatrixDef::Cmd::size) << endl;
		cout << "P2:\n" << P2->to_string(MatrixDef::Cmd::detail) << endl;

		Matrix<Pippo> P2n, P2o;
		P2o = P2n = (*P2);
		Matrix<Pippo> P2t = !(*P2);
		cout << "P2n:\n" << P2n << endl;
		cout << "P2o:\n" << P2n << endl;
		cout << "P2t:\n" << P2t << endl;

		cout << "!P2:\n" << (!(*P2)) << endl;

		Matrix<float> AR[3] = { Matrix<float>(2,2,0.1f), Matrix<float>(2,2,0.2f), Matrix<float>(2,2,0.3f) };
		for (int i = 0; i < sizeof(AR) / sizeof(AR[0]); i++)
			{
			cout << "AR[" << i << "]:\n" << AR[i] << endl;
			}

		Matrix<float> F1(2, 2, 0.5);
		cout << "F1:\n" << F1 << endl;
		Matrix<float> F2(F1);
		cout << "F2:\n" << F2 << endl;

		Matrix<Pippo> P3(*P2);
		cout << "P3:\n" << P3 << endl;

		P1.clear();
		cout << "P1:\n" << P1 << endl;

		Matrix<Pippo> P4(2, 2, Pippo(3));
		cout << "P4:\n" << P4 << endl;
		P4.clear();
		cout << "P4:\n" << P4 << endl;

		Matrix<int> J1(2, 2, 100);
		cout << "J1:\n" << J1.to_string(',', ';') << endl;

		J1.dim(4, 3);

		Matrix<int> J2(J1);
		cout << "[J2:]\t" << J2.to_string(MatrixDef::Cmd::detail) << endl;

		J1(2, 1) = 3;
		J1(3, 2) = -99;
		J1.dim(10, 10);
		cout << "J1:\n" << J1 << endl;
		J1.dim(J2);
		cout << "J1:\n" << J1 << endl;
		J1.trim(1, 2, 1, 2);
		cout << "J1:\n" << J1 << endl;
		J1.clear();

		J1.dim(5, 5);
		int i, j;
		for (i = 0; i < J1.rows(); i++)
			for (j = 0; j < J1.cols(); j++)
				J1(i, j) = rand() % 1000;
		cout << "J1:\n" << J1 << endl;

		i = 4;
		cout << "J1row[" << i << "]:\n" << (J1.get_row(i)) << endl;
		i = 3;
		cout << "J1col[" << i << "]:\n" << (J1.get_col(i)) << endl;
		cout << "J1sub:\n" << (J1.get_sub(0, 2, 0, 2)) << endl;

		J1.rem_row_col(1, 3);
		cout << "J1:\n" << J1 << endl;

		J1.rem_row(2);
		cout << "J1:\n" << J1 << endl;

		J1.transpose();
		cout << "J1tr:\n" << J1 << endl;

		J1.rem_col(1);
		cout << "J1:\n" << J1 << endl;

		#if false
		x = (*P2)(1, 3);	// Genera eccezione out of bound
		#endif

		delete P2;


		Matrix<int> m(3, 4);
		for (i = 0; i < m.rows(); i++)
			for (j = 0; j < m.cols(); j++)
				m(i, j) = rand() % 1000;
		pair<int, int> imin = m.min_row_col();
		cout << "m=\n" << m << endl << "min= " << m.min_value() << " @ " << '[' << imin.first << ',' << imin.second << ']' << endl;

		#if false
		imin = P3.min();		// Vincolo non soddisfatto
		#endif

		Matrix<int> s1(3, 4, &ramdomint);
		Matrix<int> s2(s1.rows(), s1.cols()/* se +1: errore size mismatch*/, &ramdomint);

		cout << "\ns1=\n" << s1 << endl;
		cout << "\ns2=\n" << s2 << endl;
		cout << "\ns1+s2=\n" << (s1 + s2) << endl;
		cout << "\ns1-s2=\n" << (s1 - s2) << endl;

		s1 += s2;
		cout << "\ns1+=s2;\ns1=\n" << s1 << endl;
		s1 -= s2;
		cout << "\ns1-=s2\n;s1=\n" << s1 << endl;

		Matrix<int> pr1(2, 3, ramdomint10);
		Matrix<int> pr2(3, 2, ramdomint10);
		Matrix<int> pr3 = pr1 * pr2;
		Matrix<int> pr4 = pr1 * 2;
		int ii = 2;
		Matrix<int> pr5 = ii * pr1;
		Matrix<int> pr6 = 2 * pr1;

		Matrix<long double> ld1(2, 2, 1.0001l);
		cout << "[ld1:]\t" << ld1.to_string(MatrixDef::Cmd::detail) << endl;

		Matrix<complex<double>> cp0(2, 2, cpxrc);

		Matrix<float> a(1, 3, ramdomfloat10), bb(1, 3, ramdomfloat10);



		{
		MatrixIterator<int> it(pr1);
		cout << "\npr1=\n" << pr1.to_string(MatrixDef::Cmd::detail) << endl;
		for (int *x = it.begin(); x != nullptr; x = it.next())
			{
			cout << *x << " " << *it.peek() << '\t';
			}
		}
		// Float 80 bit non supportato: long double => double in MS C++
		// __float128 bit forse supportato


		cout << "\npr2=\n" << pr2 << endl;
		cout << "\npr3=\n" << pr3 << endl;
		cout << "\npr1*pr2=\n" << (pr1 * pr2) << endl;
		cout << "\npr4=\n" << pr4 << endl;
		cout << "\npr5=\n" << pr5 << endl;
		cout << "\npr6=\n" << pr6 << endl;
		cout << "\n2*pr1=\n" << (2 * pr1) << endl;
		cout << "\npr1*2=\n" << (pr1 * 2) << endl;
		cout << "\npr6/2=\n" << (pr6 / 2) << endl;

		cout << "sizeof(char)=" << sizeof(char) << endl;
		cout << "sizeof(short)=" << sizeof(short) << endl;
		cout << "sizeof(int)=" << sizeof(int) << endl;
		cout << "sizeof(long int)=" << sizeof(long int) << endl;
		cout << "sizeof(float)=" << sizeof(float) << endl;
		cout << "sizeof(double)=" << sizeof(double) << endl;
		cout << "sizeof(long double)=" << sizeof(long double) << endl;

		cout << "\nld1=\n" << ld1 << endl;

		cout << "\ncp0=\n" << cp0 << endl;

		cout << "\na=\n" << a << endl;
		cout << "\nb=\n" << bb << endl;
		cout << "\na^bb=\n" << (a ^ bb) << endl;

		Matrix<complex<double>> idm;
		idm = Matrix<complex<double>>::Id(3);

		cout << "\nId\n" << idm << endl;

		Matrix<float> mset;
		float arrcost[3][2] = { {1.1f,2.0f},{3.0f,4.0f},{-1.0f,-2.0f} };
		mset.set(3, 2, *arrcost);
		cout << "\nmset.to_string:\n" << mset << endl;

		Matrix<float> mset1;
		// mset1.set(2, 2,{{10,20},{-3,-4}});
		mset1.set(2, 2, {	10.0f,	20.0f,
							-3.0f,	-4.0f
						});
		cout << "\nmset1.to_string:\n" << mset1 << endl;

		// mset1.set(2, 2, { 10,20,-3,-4,5});	// Error, sizes do not match

		LinearSys<float, float> *ls1 = new LinearSys<float, float>();
		LinearSys<double, double> *ls2 = new LinearSys<double, double>();
		LinearSys<complex<float>, float> *ls3 = new LinearSys<complex<float>, float>();
		LinearSys<complex<double>, double> *ls4 = new LinearSys<complex<double>, double>();
		LinearSys<complex<long double>, long double> *ls5 = new LinearSys<complex<long double>, long double>();
		ls4->set_eps_zero_ratio(10);
		cout << "\nls1->to_string:\n" << ls1->to_string() << endl;
		cout << "\nls2->to_string:\n" << ls2->to_string() << endl;
		cout << "\nls3->to_string:\n" << ls3->to_string() << endl;
		cout << "\nls4->to_string:\n" << ls4->to_string() << endl;
		cout << "\nls5->to_string:\n" << ls5->to_string() << endl;

		// LinearSys<float,string> *ls2 = new LinearSys<float, string>();	// Error C7500: nessuna funzione soddisfa i vincoli

		//Matrix<double> A(3,3,ramdomdouble20);

		Matrix<double> A, xo, b, x;


		A.set(3, 3, {	1, 2, 3,
						-2, 4, -7,
						5, 6, 1 });
		xo.set(3, 1, {	1, -5, 8 });


		/*
		A.set(4, 4, { 0, 2, 0, -1,
									2, -1, 1, -2,
									1, 0, -2, 1,
									-1, 3, 1, 1
									});
		xo.set(4, 1, {1, 0, -3, 3});
		*/

		b = A * xo;

		cout << "\nA:\n" << A << endl;
		cout << "\nxo:\n" << xo << endl;
		cout << "\nb=A*xo:\n" << b << endl;

		ls2->factor(A);

		cout << "\nls2:Factor()\n" << ls2->to_string() << endl;

		/*	Verifica con GNU Octave
		>> A = [1 2 3; -2 4 -7; 5 6 1]
		A =		1   2   3
				-2   4  -7
				5   6   1
		>> [L,U]=lu(A)
		L =		0.2000	0.1250	1.0000
				-0.4000	1.0000	0
				1.0000	0		0
		U =		5.0000	6.0000	1.0000
				0	6.4000	-6.6000
				0	0	3.6250
		>> bb=[5 3 10]'
		bb =		5
				3
				10
		>> x=A\bb
		x =		-0.1552
				1.6983
				0.5862
		>> A*x
		ans =	5
				3
				10
		*/

		if (ls2->solve_check())
			{
			ls2->solve(x, b);
			} else
			{
			cout << "\nsolve_check non superato" << endl;
			}
			cout << "\nx=\n" << x << endl;
			cout << "\nA*x=\n" << (A * x) << endl;
			cout << "\nx-xo\n" << (x - xo) << endl;


		}
	catch (const std::runtime_error ex)
		{
		cout << ex.what();
		stop_msg();
		exit(1);
		}

		stop_msg();
	}

void stop_msg()
	{
	cout << endl << "<enter> per chiudere" << endl;
	char tempchar = getchar();
	}

static int ramdomint(int i, int j)
	{
	return 1 + rand() % 100;
	}

static int ramdomint10(int i, int j)
	{
	return (rand() % 10) - (rand() % 5);
	}

static float ramdomfloat10(int i, int j)
	{
	return (float)((rand() % 10) - (rand() % 5));
	}

static double ramdomdouble20(int i, int j)
	{
	return (double)((rand() % 20) - (rand() % 10));
	}

static std::complex<double> cpxrc(int r, int c)
	{
	return std::complex<double>(r * 10 + ramdomint10(0, 0), c * 10 + ramdomint10(0, 0));
	}

