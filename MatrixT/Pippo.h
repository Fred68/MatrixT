
#ifndef PIPPO
#define PIPPO

#define ERR_NO_ASSIGN			// Prova vincolo operatore assegnazione
#undef ERR_NO_ASSIGN

#include <iostream>
using std::cout;
using std::endl;
using std::ostream;

class Pippo
	{
	protected:
		int _p = -1;
		double _p2;
		int _in[3] = {1,2,3};
	public:
		Pippo()
			{
			_p = 0;
			_p2 = (double)_p / 2.0;
			#ifdef _DEBUG
			cout << "ctor Pippo()" << endl;
			#endif
			}
		Pippo(int x)
			{
			_p = x;
			#ifdef _DEBUG
			cout << "ctor Pippo(int x)" << endl;
			#endif
			}
		int Get() { return _p; }
		void Set(int p) { _p = p; }
		friend ostream &operator<<(ostream &os, const Pippo &pp);
		Pippo &operator=(const Pippo &pp)
		#ifdef ERR_NO_ASSIGN
			= delete;
		#else
			{
			if (this != &pp)
				{
				_p = pp._p;
				#ifdef _DEBUG
				cout << "Pippo operator=()" << endl;
				#endif
				}
			else
				{
				#ifdef _DEBUG
				cout << "Pippo auto assignment in operator=()" << endl;
				#endif
				}
			return *this;
			}
		#endif
	};

#endif

