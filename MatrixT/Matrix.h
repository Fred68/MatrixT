

#ifndef MATRIXH
#define MATRIXH

/*************************************/
/*                                   */
/*  Matrix.h                         */
/*                                   */
/*************************************/


#include <iostream>
#include <concepts>
#include <utility>
#include <string>
#include <sstream>
#include <complex>
#include <limits>
#include <vector>
#include <cmath>
#include <cstddef>

// #define SPECIALITERATOR

// #include <type_traits>	// Per complex
using std::swap;
using std::string;
using std::runtime_error;
using std::pair;
using std::is_floating_point;
using std::numeric_limits;
using std::vector;
using std::stringstream;
using std::endl;
using std::ostream;

#ifdef _DEBUG
using std::cout;
#endif




namespace matrix
	{

	#pragma region PROTOTYPING

	class MatrixDef;
	template<typename DATA> class Matrix;
	template<typename DATA> class MatrixIterator;	/* Prototyping della classe generica che deve essere friend */

	template<typename U> ostream& operator<<(ostream& ostr, Matrix<U>& m);
	template<typename U> ostream& operator<<(ostream& ostr, const Matrix<U>& m);

	#pragma endregion

	#pragma region CONCEPTS

	/* Assegnazione semplice */
	template<typename T>
	concept RQassign = requires (T x, T y) { x = y; };

	/* Assegnazione e ordinamento*/
	template<typename T>
	concept RQorder = RQassign<T> && std::totally_ordered<T>;

	/* Somma (valida anche per stringhe) */
	template<typename T>
	concept RQsum = RQassign<T> && requires (T x, T y) { x += y; x + y; };

	/* Somma, differerenza, prodotto */
	template<typename T>
	concept RQsumdifprod = RQsum<T> && requires (T x, T y) { x -= y; x - y; x *= y; x * y; };

	/* Per tutti i numeri a virgola mobile o complessi */
	template<typename T>
	concept RQfloat = RQsumdifprod<T> && (std::is_floating_point<T>::value || std::is_same_v<T, std::complex<typename T::value_type>>);

	/* Per tutti i numeri a virgola mobile o complessi con funzione abs()*/
	template<typename T, class A>
	concept RQfloatabs = RQfloat<T> && RQfloat<A> && requires (T x, A a) { a = abs(x); };

	#pragma endregion


	class MatrixDef
	{
		public:
			enum class Cmd { size, detail };
			inline static const string ERR_ALLOC = "Allocation failed";
			inline static const string ERR_OUTOFBOUND = "Subscript out of bounds";
			inline static const string ERR_WRONGPARAM = "Wrong parameter";
			inline static const string ERR_SZMISMATCH = "Sizes do not match";
			inline static const string ERR_NOTVECTOR = "Not a vector";
			inline static const string ERR_ZEROSIZE = "Zero size";
			inline static const string ERR_NOTSQUARE = "Not square";
			inline static const string ERR_SINGULAR = "Singular";
			inline static const string ERR_PIVOT = "Wrong pivot";
	};


	#ifdef SPECIALITERATOR
	/********************************************************/
	// PARTE RIFATTA CON SEMPLICE ITERATORE STANDARD
	/********************************************************/
	template<typename DATA> class MatrixIterator
		{
		private:
			Matrix<DATA> *_m;
			int ir, ic;
		public:
			/* Ctor con argomenti */
			MatrixIterator(Matrix<DATA> &m)
				{
				_m = &m;
				ir = ic = 0;
				_m->_iterators++;
				};
			/* Dtor con argomenti */
			~MatrixIterator()
				{
				_m->_iterators--;
				}
			/* Iteration */
			DATA *begin()
				{
				ir = ic = 0;
				if ((_m->_row == 0) || (_m->_col == 0))
					{
					throw std::runtime_error(MatrixDef::ERR_ZEROSIZE);
					}
				return _m->dat;
				}
			DATA *end()
				{
				ir = _m->_row - 1;
				ic = _m->_col - 1;
				if ((ir < 0) || (ic < 0))
					{
					throw std::runtime_error(MatrixDef::ERR_ZEROSIZE);
					}
				return _m->dat[ir * _m->_col + ic];
				}
			DATA *next()
				{
				ic++;
				if (ic >= _m->_col)
					{
					ic = 0;
					ir++;
					}
				if (ir >= _m->_row)				// Oltre fine matrice
					{
					ic = ir = 0;				// Azzera l'iteratore
					return (DATA*) nullptr;
					}
				return _m->dat + ir * _m->_col + ic;
				}
			DATA *peek()
				{
				return _m->dat + ir * _m->_col + ic;
				}
		};
	#endif	
	/********************************************************/

	template<typename DATA> class Matrix
		{
		private:
			int _row;							// Numeri di righe...
			int _col;							// ...e colonne
			DATA *dat;							// Puntatore ai dati
			static DATA *_empty;				// Dato vuoto
			#ifdef SPECIALITERATOR
			int _iterators;						// Numero di iteratori
			#endif
			static size_t _datasize;			// DATA size
		public:
			struct iterator
			{
				using iterator_category = std::forward_iterator_tag;	// ->> Da rivedere con sintassi C++20
				explicit iterator(DATA* ptr) : _ptr{ ptr } {}
				DATA& operator*() const { return *_ptr; }
				DATA* operator->() const { return _ptr; }
				iterator& operator++() { _ptr++; return *this; }
				iterator operator++(int) { iterator tmp = *this; _ptr++; return tmp; }
				friend bool operator!= (const iterator& a, const iterator& b) { return a._ptr != b._ptr; };
			private:
				DATA* _ptr;
			};

			/* Ctor */
			Matrix()
			{
				_row = _col = 0;
				dat = (DATA*)nullptr;
				#ifdef SPECIALITERATOR
				_iterators = 0;
				#endif
#ifdef _DEBUG
				cout << "Matrix()" << endl;
#endif
			}
			Matrix(int rows, int cols)
			{
				if ((rows > 0) && (cols > 0))
				{
					_row = rows;
					_col = cols;
					dat = new DATA[rows * cols];
					if (!dat)							// Se fallita allocazione, esce con errore
					{
						_row = _col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
				}
				else
				{									// Matrice vuota
					_row = _col = 0;
					dat = (DATA*)nullptr;
				}
				#ifdef SPECIALITERATOR
				_iterators = 0;
				#endif
#ifdef _DEBUG
				cout << "Matrix(int rows, int cols)" << endl;
#endif
			}
			Matrix(int rows, int cols, DATA d) requires RQassign<DATA>
			{
#pragma warning (disable : 6385)
				if ((rows > 0) && (cols > 0))
				{
					_row = rows;
					_col = cols;
					dat = new DATA[rows * cols];
					if (!dat)
					{
						_row = _col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
					else
					{
						int ir, ic;
						for (ir = 0; ir < _row; ir++)
							for (ic = 0; ic < _col; ic++)
								dat[ir * _col + ic] = d;
					}
				}
				else
				{									// Matrice vuota
					_row = _col = 0;
					dat = (DATA*)nullptr;
				}
				#ifdef SPECIALITERATOR
				_iterators = 0;
				#endif
#ifdef _DEBUG
				cout << "Matrix(int rows, int cols, DATA d)" << endl;
#endif
#pragma warning (default : 6385)
			}
			Matrix(int rows, int cols, DATA(*pf)(int r, int c)) requires RQassign<DATA>
			{
				if ((rows > 0) && (cols > 0))
				{
					_row = rows;
					_col = cols;
					dat = new DATA[rows * cols];
					if (!dat)
					{
						_row = _col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
					else
					{
						int ir, ic;
						for (ir = 0; ir < _row; ir++)
							for (ic = 0; ic < _col; ic++)
								dat[ir * _col + ic] = (*pf)(ir, ic);
					}
				}
				else
				{									// Matrice vuota
					_row = _col = 0;
					dat = (DATA*)nullptr;
				}
#ifdef SPECIALITERATOR
				_iterators = 0;
#endif
#ifdef _DEBUG
				cout << "Matrix(int rows, int cols, DATA d)" << endl;
#endif
			}
			/* Copy & move ctor */
			Matrix(const Matrix& m) requires RQassign<DATA> : _row{m._row}, _col{m._col}
			#ifdef SPECIALITERATOR
			,_iterators{0}
			#endif
			{
				if ((_row > 0) && (_col > 0))
				{
					dat = new DATA[_row * _col];
					if (dat)
					{
						int ir, ic;							// Non usa memcpy((void*)dat, (void*)m.dat, row * col * sizeof(DATA)); ma '='.
						for (ir = 0; ir < _row; ir++)		// Ricopia i valori. Usa operatore di assegnazione.
							for (ic = 0; ic < _col; ic++)
								dat[ir * _col + ic] = m.dat[ir * _col + ic];
					}
					else
					{
						_row = _col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
				}
				else
				{
					_row = _col = 0;
					dat = (DATA*)nullptr;
				}
				
#ifdef _DEBUG
				cout << "Matrix(const Matrix& m)" << endl;
#endif
			}
			Matrix(Matrix&& m) : _row{m._row}, _col{m._col}, dat{m.dat}
#ifdef SPECIALITERATOR
			, _iterators{ 0 }
#endif

			{
				m.dat = nullptr;
				m._row = m._col = 0;
				m._iterators = 0;
#ifdef _DEBUG
				cout << "Matrix(Matrix&& m)" << endl;
#endif
			}
			// Copy and move assignment
			Matrix<DATA> &operator=(const Matrix<DATA> &m) requires RQassign<DATA>
			{
				if (this != &m)
				{
					DATA *tmp;
					_row = m._row;
					_col = m._col;
					if ((_row > 0) && (_col > 0))
					{
						tmp = new DATA[_row * _col];
						if (tmp)
						{
							int ir, ic;							// Non usa memcpy((void*)dat, (void*)m.dat, row * col * sizeof(DATA)); ma '='.
							for (ir = 0; ir < _row; ir++)		// Ricopia i valori. Usa operatore di assegnazione.
								for (ic = 0; ic < _col; ic++)
									tmp[ir * _col + ic] = m.dat[ir * _col + ic];
							delete[] dat;
							dat = tmp;
						}
						else
						{
							_row = _col = 0;
							throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					}
					else
					{									// Matrice vuota
						_row = _col = 0;
						dat = (DATA*)nullptr;
					}
				}
#ifdef _DEBUG
				cout << "Matrix::&operator=(const Matrix&)" << endl;
#endif
				return *this;
			}
			Matrix<DATA> &operator=(Matrix<DATA> &&m)
			{
				if(this != &m)
				{
					_col = m._col;
					_row = m._row;
					#ifdef SPECIALITERATOR
					_iterators = m-_iterators;
					#endif
					dat = m.dat;
					m.dat = nullptr;
					m._row = m._col = 0;
					m._iterators = 0;
				}
#ifdef _DEBUG
				cout << "Matrix::&operator=(Matrix&&)" << endl;
#endif
				return *this;
			}
			/* Dtor */
			~Matrix()
			{
				if (dat != (DATA*)nullptr)
				{
					delete[] dat;
				}
#ifdef SPECIALITERATOR
				if (_iterators > 0)
				{
					// throw std::runtime_error(MatrixDef::ERR_ALLOC);	// Non interrompere un distruttore.
#ifdef _DEBUG
					cout << "~Matrix() with pending " << _iterators << " iterators." << endl;
#endif
				}
#endif // SPECIALITERATOR
#ifdef _DEBUG
				cout << "~Matrix()" << endl;
#endif
			}
			// Accesso in lettura
			int rows()
				{
				return _row;
				}
			int cols()
				{
				return _col;
				}
			// Accesso ai dati
			DATA &operator()(int irow, int icol)
				{
				if ((irow < _row) && (irow >= 0) && (icol < _col) && (icol >= 0))
					{
					return dat[irow * _col + icol];
					}
				else
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				}
			Matrix<DATA> &get_row(int irow) requires RQassign<DATA>
				{
				int ic;
				if ((irow < 0) || (irow >= _row))
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				Matrix<DATA> *m = new Matrix<DATA>(1, _col);
				for (ic = 0; ic < _col; ic++)
					{
					m->dat[ic] = dat[irow * _col + ic];
					}
				#ifdef _DEBUG
				cout << "Matrix::get_row(int irow)" << endl;
				#endif
				return *m;
				}
			Matrix<DATA> &get_col(int icol) requires RQassign<DATA>
				{
				int ir;
				if ((icol < 0) || (icol >= _col))
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				Matrix<DATA> *m = new Matrix<DATA>(_row, 1);
				for (ir = 0; ir < _row; ir++)
					{
					m->dat[ir] = dat[ir * _col + icol];
					}
				#ifdef _DEBUG
				cout << "Matrix::get_col(int icol)" << endl;
				#endif
				return *m;
				}
			Matrix<DATA> &get_sub(int row_ini, int col_ini, int n_row, int n_col) requires RQassign<DATA>
				{
				if ((n_row < 0) || (n_col < 0))
					{
					throw std::runtime_error(MatrixDef::ERR_WRONGPARAM);
					}
				int row_fin = row_ini + n_row - 1;
				int col_fin = col_ini + n_col - 1;
				if ((row_ini < 0) || (row_fin >= _row) || (col_ini < 0) || (col_fin >= _col))
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				Matrix<DATA> *m = new Matrix<DATA>(row_fin - row_ini + 1, col_fin - col_ini + 1);
				if ((m->_row > 0) && (m->_col > 0))
					{
					int ir, ic;
					for (ir = 0; ir < n_row; ir++)
						{
						for (ic = 0; ic < n_col; ic++)
							{
							m->dat[ir * m->_col + ic] = dat[(row_ini + ir) * _col + (col_ini + ic)];
							}
						}
					}
				#ifdef _DEBUG
				cout << "Matrix::get_sub(int row_ini, int col_ini, int n_row, int n_col)" << endl;
				#endif
				return *m;
				}
			// set
			void clear()
				{
				if (dat != (DATA*)nullptr)
					{
					delete[] dat;
					dat = (DATA*)nullptr;
					}
				_col = _row = 0;
				#ifdef _DEBUG
				cout << "Clear()" << endl;
				#endif
				}
			void set(int rows, int cols, DATA *d) requires RQassign<DATA>
				{
				DATA *datnew;							// Puntatore ai nuovi dati
				if ((rows < 0) || (cols < 0))			// Verifica le nuove dimensioni
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				if ((rows == 0) || (cols == 0))			// Se una dimensione è nulla
					{
					if (dat) delete[] dat;				// Dealloca la vecchia matrice
					_row = _col = 0;					// e azzera
					dat = (DATA*)nullptr;
					return;
					}
				// if (array.size() != rows * cols) {throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);}
				datnew = new DATA[rows * cols];			// Alloca nuova matrice
				if (!datnew)							// Verifica allocazione avvenuta
					{
					throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
				int ir, ic;								// Ricopia i valori dal vettore alla matrice
				for (ir = 0; ir < rows; ir++)
					for (ic = 0; ic < cols; ic++)
						datnew[ir * cols + ic] = d[ir * cols + ic];
				if (dat) delete[] dat;					// Dealloca la vecchia matrice
				_row = rows;							// Imposta i nuovi indici
				_col = cols;
				dat = datnew;							// Imposta il nuovo puntatore

				#ifdef _DEBUG
				cout << "void Matrix::set(int rows, int cols, DATA *d). ---> No array boundary check !" << endl;
				#endif
				return;
				}
			void set(int rows, int cols, const std::vector<DATA>& array) requires RQassign<DATA>
				{
				DATA *datnew;							// Puntatore ai nuovi dati
				if ((rows < 0) || (cols < 0))			// Verifica le nuove dimensioni
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				if ((rows == 0) || (cols == 0))			// Se una dimensione è nulla
					{
					if (dat) delete[] dat;				// Dealloca la vecchia matrice
					_row = _col = 0;					// e azzera
					dat = (DATA*)nullptr;
					return;
					}
				if (array.size() != rows * cols)
					{
					throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
					}
				datnew = new DATA[rows * cols];			// Alloca nuova matrice
				if (!datnew)							// Verifica allocazione avvenuta
					{
					throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
				int ir, ic;								// Ricopia i valori dal vettore alla matrice
				for (ir = 0; ir < rows; ir++)
					for (ic = 0; ic < cols; ic++)
						datnew[ir * cols + ic] = array[ir * cols + ic];
				if (dat) delete[] dat;					// Dealloca la vecchia matrice
				_row = rows;							// Imposta i nuovi indici
				_col = cols;
				dat = datnew;							// Imposta il nuovo puntatore

				#ifdef _DEBUG
				cout << "void Matrix::set(int rows, int cols, const std::vector<DATA>& array)" << endl;
				#endif
				return;
				}
			// Trasposta
			Matrix<DATA> &operator!() requires RQassign<DATA>
				{
				Matrix<DATA> *tmp = new Matrix<DATA>();							// Alloca nuovo Matrix vuoto
				int ir, ic;
				if ((_row > 0) && (_col > 0))
					{
					tmp->dat = new DATA[_row * _col];							// Alloca spazio.
					if (tmp->dat)
						{
						for (ir = 0; ir < _row; ir++)							// Ricopia i valori, scambiandone le posizioni
							for (ic = 0; ic < _col; ic++)
								tmp->dat[ic * _row + ir] = dat[ir * _col + ic];
						tmp->_row = _col;
						tmp->_col = _row;
						} else
						{
						tmp->_row = tmp->_col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					} else
					{																	// Matrice nulla
					tmp->_row = tmp->_col = 0;
					tmp->dat = (DATA*)nullptr;
					}
					#ifdef _DEBUG
					cout << "Matrix::&operator!()" << endl;
					#endif
					return *tmp;
				}
			void transpose() requires RQassign<DATA>				// Traspone
				{
				DATA *tmp = new DATA[_row * _col];									// Alloca spazio.
				int ir, ic;
				if (!tmp)
					{
					throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
				for (ir = 0; ir < _row; ir++)										// Ricopia i valori, scambiandone le posizioni
					for (ic = 0; ic < _col; ic++)
						tmp[ic * _row + ir] = dat[ir * _col + ic];
				ic = _row;															// Scambia righe e colonne	
				_row = _col;
				_col = ic;
				delete[] dat;														// Dealloca i dati vecchi
				dat = tmp;															// Riassegna il puntatore ai dati nuovi
				#ifdef _DEBUG
				cout << "Matrix::transpose()" << endl;
				#endif
				return;
				}
			// To String() 	Operazioni su stream generano link error con: friend ostream &operator<<(ostream &stream, const Matrix &m);	
			string to_string(char col_sep = '\t', char row_sep = '\n')
				{
				int ir, ic;
				std::stringstream ss;
				ss << "[";
				for (ir = 0; ir < _row; ir++)
					{
					for (ic = 0; ic < _col; ic++)
						{
						ss << dat[ir * _col + ic];
						if (ic != _col - 1)	ss << col_sep;
						}
					ss << row_sep;
					}
				ss << "] R" << _row << " x C" << _col;
				#ifdef _DEBUG
				cout << "Matrix::ToString(char col_sep, char row_sep)" << endl;
				#endif
				return ss.str();
				}
			string to_string(MatrixDef::Cmd cmd)
				{
				std::stringstream ss;
				switch (cmd)
					{
					case MatrixDef::Cmd::size:
					ss << "R" << _row << " x C" << _col;
					break;
					case MatrixDef::Cmd::detail:
					ss << "R" << _row << "xC" << _col;
#ifdef SPECIALITERATOR
					ss << ",iters=" << _iterators;
#endif
					ss << ",dtsz=" << _datasize << ",empty=" << *_empty;
					break;
					default:
					break;
					}

				#ifdef _DEBUG
				cout << "Matrix::ToString(int cmd)" << endl;
				#endif
				return ss.str();
				}
			// Modifica dimensioni
			bool dim(int rows, int cols) requires RQassign<DATA>
				{
				DATA *datnew;							// Puntatore ai nuovi dati
				int ir, ic;
				if ((rows < 0) || (cols < 0))			// Verifica le nuove dimensioni
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				if ((rows == 0) || (cols == 0))			// Se una dimensione è nulla
					{
					if (dat) delete[] dat;				// Delalloca la vecchia matrice
					_row = _col = 0;					// e azzera
					dat = (DATA*)nullptr;
					return false;
					}
				datnew = new DATA[rows * cols];			// Alloca nuova matrice
				if (!datnew)							// Verifica allocazione avvenuta
					{
					throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
				for (ir = 0; ir < rows; ir++)			// Ricopia la vecchia matrice nella nuova
					for (ic = 0; ic < cols; ic++)
						{
						if ((ir < _row) && (ic < _col))
							datnew[ir * cols + ic] = dat[ir * _col + ic];	// Ricopia valori...
						else
							datnew[ir * cols + ic] = (DATA)0;		// oppure azzera se fuori indice
						}
				if (dat) delete[] dat;					// Dealloca la vecchia matrice
				_row = rows;							// Imposta i nuovi indici
				_col = cols;
				dat = datnew;							// Imposta il nuovo puntatore
				#ifdef _DEBUG
				cout << "Matrix::dim(int rows, int cols)" << endl;
				#endif
				return true;
				}
			bool dim(Matrix<DATA> &m)  requires RQassign<DATA>
				{
				return dim(m.rows(), m.cols());
				}
			bool trim(int row1, int row2, int col1, int col2) requires RQassign<DATA>
				{
				// Mantiene indici tra row1,row2 e col1,col2 compresi.
				// Se row1>row2 o col1>col2, azzera la matrice;
				DATA *datnew;							// Puntatore a nuovi dati
				int rownew, colnew;
				int ir, ic;

				rownew = row2 - row1 + 1;					// Calcola nuove dimensioni
				colnew = col2 - col1 + 1;
				if ((colnew < 0) || (rownew < 0))			// Indici errati
					{
					throw std::runtime_error(MatrixDef::ERR_WRONGPARAM);
					}

				if ((colnew <= 0) || (rownew <= 0))			// Se nulle, azzera la matrice
					{
					if (dat) delete[] dat;
					_row = _col = 0;
					dat = (DATA*)nullptr;
					return false;
					}
				if ((row1 < 0) || (row2 < 0) || (row1 >= _row) || (row2 >= _row))	// Verifica correttezza indici
					return false;
				if ((col1 < 0) || (col2 < 0) || (col1 >= _col) || (col2 >= _col))
					return false;
				if (row1 > row2)							// Ordina gli indici
					ir = row1, row1 = row2, row2 = ir;
				if (col1 > col2)
					ic = col1, col1 = col2, col2 = ic;

				datnew = new DATA[(rownew) * (colnew)];	// Alloca nuova matrice
				if (!datnew)							// Verifica avvenuta allocazione
					{
					throw std::runtime_error(MatrixDef::ERR_ALLOC);
					}
				for (ir = row1; ir <= row2; ir++)				// Ricopia i valori
					for (ic = col1; ic <= col2; ic++)
						{
						datnew[(ir - row1) * colnew + (ic - col1)] = dat[ir * _col + ic];
						}
				if (dat) delete[] dat;					// Dealloca la vecchia matrice
				_row = rownew;							// Imposta i nuovi indici
				_col = colnew;
				dat = datnew;							// Imposta il nuovo puntatore
				#ifdef _DEBUG
				cout << "Matrix::trim(int row1, int row2, int col1, int col2)" << endl;
				#endif
				return true;
				}
			bool rem_row_col(int irow, int icol) requires RQassign<DATA>
				{
				DATA *datnew;							// Puntatore a nuovi dati
				int rownew, colnew;						// Nuove dimensioni
				int ir, ic;								// Indici nuova matrice
				int irold, icold;						// Indici vecchia matrice

				if ((irow < 0) || (irow >= _row) || (icol < 0) || (icol >= _col))		// Verifica correttezza indici
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				rownew = _row - 1;						// Nuove dimensioni
				colnew = _col - 1;

				if ((colnew != 0) && (rownew != 0))		// Se matrice di dimensione non nulla
					{
					datnew = new DATA[rownew * colnew];	// Alloca nuova matrice (row-1, col-1)
					} else
					{
					if (dat)		delete[] dat;			// Se dimensione nulla, dealloca la vecchia
					dat = (DATA*)nullptr;				// Azzera la matrice
					_col = _row = 0;
					return false;
					}
					if (!datnew)							// Se fallita allocazione
						{
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					for (ir = 0, irold = 0; ir < rownew; ir++, irold++)
						{
						if (irold == irow)	irold++;		// Salta all'indice successivo, se uguale alla riga da elimonare
						for (ic = 0, icold = 0; ic < colnew; ic++, icold++)
							{
							if (icold == icol)	icold++;	// Salta all'indice successivo, se uguale alla colonna da elimonare
							datnew[ir * colnew + ic] = dat[irold * _col + icold];	// Copia il valore
							}
						}
					if (dat)		delete[] dat;
					_row = rownew;
					_col = colnew;
					dat = datnew;
					#ifdef _DEBUG
					cout << "Matrix::rem_row_col(int irow, int icol)" << endl;
					#endif
					return true;
				}
			bool rem_row(int irow) requires RQassign<DATA>
				{
				DATA *datnew;							// Puntatore a nuovi dati
				int rownew, colnew;						// Nuove dimensioni
				int ir, ic;								// Indici nuova matrice
				int irold, icold;						// Indici vecchia matrice

				if ((irow < 0) || (irow >= _row))		// Verifica correttezza indici
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				rownew = _row - 1;						// Nuove dimensioni
				colnew = _col;

				if ((colnew != 0) && (rownew != 0))		// Se matrice di dimensione non nulla
					{
					datnew = new DATA[rownew * colnew];	// Alloca nuova matrice (row-1, col-1)
					} else
					{
					if (dat)		delete[] dat;		// Se dimensione nulla, dealloca la vecchia
					dat = (DATA*)nullptr;				// Azzera la matrice
					_col = _row = 0;
					return false;
					}
					if (!datnew)							// Se fallita allocazione
						{
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					for (ir = 0, irold = 0; ir < rownew; ir++, irold++)
						{
						if (irold == irow)	irold++;		// Salta all'indice successivo, se uguale alla riga da elimonare
						for (ic = 0, icold = 0; ic < colnew; ic++, icold++)
							{
							datnew[ir * colnew + ic] = dat[irold * _col + icold];	// Copia il valore
							}
						}
					if (dat)		delete[] dat;
					_row = rownew;
					_col = colnew;
					dat = datnew;
					#ifdef _DEBUG
					cout << "Matrix::rem_row(int irow)" << endl;
					#endif
					return true;
				}
			bool rem_col(int icol) requires RQassign<DATA>
				{
				DATA *datnew;							// Puntatore a nuovi dati
				int rownew, colnew;						// Nuove dimensioni
				int ir, ic;								// Indici nuova matrice
				int irold, icold;						// Indici vecchia matrice

				if ((icol < 0) || (icol >= _col))		// Verifica correttezza indici
					{
					throw std::runtime_error(MatrixDef::ERR_OUTOFBOUND);
					}
				rownew = _row;							// Nuove dimensioni
				colnew = _col - 1;

				if ((colnew != 0) && (rownew != 0))		// Se matrice di dimensione non nulla
					{
					datnew = new DATA[rownew * colnew];	// Alloca nuova matrice (row-1, col-1)
					} else
					{
					if (dat)		delete[] dat;		// Se dimensione nulla, dealloca la vecchia
					dat = (DATA*)nullptr;				// Azzera la matrice
					_col = _row = 0;
					return false;
					}
					if (!datnew)							// Se fallita allocazione
						{
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					for (ir = 0, irold = 0; ir < rownew; ir++, irold++)
						{
						for (ic = 0, icold = 0; ic < colnew; ic++, icold++)
							{
							if (icold == icol)	icold++;	// Salta all'indice successivo, se uguale alla colonna da elimonare
							datnew[ir * colnew + ic] = dat[irold * _col + icold];	// Copia il valore
							}
						}
					if (dat)	delete[] dat;
					_row = rownew;
					_col = colnew;
					dat = datnew;
					#ifdef _DEBUG
					cout << "Matrix::rem_col(int icol)" << endl;
					#endif
					return true;
				}
			// Min, max
			pair<int,int> &min_row_col() requires RQorder<DATA>
				{
				pair<int, int> *r = new pair<int, int>(-1, -1);
				DATA tmp;
				if ((_row > 0) && (_col > 0))
					{
					int ir, ic;
					tmp = dat[0];
					ir = ic = 0;
					for (ir = 0; ir < _row; ir++)
						{
						for (ic = 0; ic < _col; ic++)
							{
							if (dat[ir * _col + ic] < tmp)
								{
								tmp = dat[ir * _col + ic];
								r->first = ir;
								r->second = ic;
								}
							}
						}
					}
				#ifdef _DEBUG
				cout << "Matrix::min_row_col()" << endl;
				#endif
				return *r;
				}
			pair<int, int> &max_row_col() requires RQorder<DATA>
				{
				pair<int, int> *r = new pair<int, int>(-1, -1);
				DATA tmp;
				if ((_row > 0) && (_col > 0))
					{
					int ir, ic;
					tmp = dat[0];
					ir = ic = 0;
					for (ir = 0; ir < _row; ir++)
						{
						for (ic = 0; ic < _col; ic++)
							{
							if (dat[ir * _col + ic] > tmp)
								{
								tmp = dat[ir * _col + ic];
								r->first = ir;
								r->second = ic;
								}
							}
						}
					}
				#ifdef _DEBUG
				cout << "Matrix::max_row_col()" << endl;
				#endif
				return *r;
				}
			DATA &min_value() requires RQorder<DATA>
				{
				pair<int,int> i = min_row_col();	// Occasionale errore su funzione min() ?
				return (*this)(i.first, i.second);
				}
			DATA &max_value() requires RQorder<DATA>
				{
				pair<int, int> i = max_row_col();
				return (*this)(i.first, i.second);
				}
			// Operatore somma
			Matrix <DATA> &operator+(const Matrix <DATA> &m) requires RQsum<DATA>
				{
				if ((m._row != _row) || (m._col != _col))
					{
					throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
					}
				Matrix<DATA> *tmp = new Matrix<DATA>();							// Alloca nuova Matrix vuota
				int ir, ic;
				if ((_row > 0) && (_col > 0))
					{
					tmp->dat = new DATA[_row * _col];							// Alloca spazio.
					if (tmp->dat)
						{
						for (ir = 0; ir < _row; ir++)							// Somma
							for (ic = 0; ic < _col; ic++)
								tmp->dat[ir * _col + ic] = dat[ir * _col + ic] + m.dat[ir * _col + ic];
						tmp->_row = _row;
						tmp->_col = _col;
						} else
						{
						tmp->_row = tmp->_col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					} else
					{															// Matrice nulla
					tmp->_row = tmp->_col = 0;
					tmp->dat = (DATA*)nullptr;
					}
					#ifdef _DEBUG
					cout << "Matrix::operator+(const Matrix <DATA> &m)" << endl;
					#endif
					return *tmp;
				}
			Matrix <DATA> &operator+=(const Matrix <DATA> &m) requires RQsum<DATA>
				{
				if ((m._row != _row) || (m._col != _col))
					{
					throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
					}
				int ir, ic;														// Non alloca nulla, non deve cambiare dimensioni della matrice.
				if ((_row > 0) && (_col > 0))
					{
					for (ir = 0; ir < _row; ir++)								// Somma
						for (ic = 0; ic < _col; ic++)
							dat[ir * _col + ic] = dat[ir * _col + ic] + m.dat[ir * _col + ic];
					}
				#ifdef _DEBUG
				cout << "Matrix::operator+=(const Matrix <DATA> &m)" << endl;
				#endif
				return *this;
				}
			// Differenza e prodotto di matrici
			Matrix <DATA> &operator-(const Matrix <DATA> &m) requires RQsumdifprod<DATA>
				{
				if ((m._row != _row) || (m._col != _col))
					{
					throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
					}
				Matrix<DATA> *tmp = new Matrix<DATA>();							// Alloca nuova Matrix vuota
				int ir, ic;
				if ((_row > 0) && (_col > 0))
					{
					tmp->dat = new DATA[_row * _col];							// Alloca spazio.
					if (tmp->dat)
						{
						for (ir = 0; ir < _row; ir++)							// Somma
							for (ic = 0; ic < _col; ic++)
								tmp->dat[ir * _col + ic] = dat[ir * _col + ic] - m.dat[ir * _col + ic];
						tmp->_row = _row;
						tmp->_col = _col;
						} else
						{
						tmp->_row = tmp->_col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					} else
					{															// Matrice nulla
					tmp->_row = tmp->_col = 0;
					tmp->dat = (DATA*)nullptr;
					}
					#ifdef _DEBUG
					cout << "Matrix::operator-(const Matrix <DATA> &m)" << endl;
					#endif
					return *tmp;
				}
			Matrix <DATA> &operator-=(const Matrix <DATA> &m) requires RQsumdifprod<DATA>
				{
				if ((m._row != _row) || (m._col != _col))
					{
					throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
					}
				int ir, ic;														// Non alloca nulla, non deve cambiare dimensioni della matrice.
				if ((_row > 0) && (_col > 0))
					{
					for (ir = 0; ir < _row; ir++)								// Somma
						for (ic = 0; ic < _col; ic++)
							dat[ir * _col + ic] = dat[ir * _col + ic] - m.dat[ir * _col + ic];
					}
				#ifdef _DEBUG
				cout << "Matrix::operator-=(const Matrix <DATA> &m)" << endl;
				#endif
				return *this;
				}
			Matrix <DATA> &operator*(const Matrix <DATA> &m) requires RQsumdifprod<DATA>
				{
				if (_col != m._row)
					{
					throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
					}
				Matrix<DATA> *tmp = new Matrix<DATA>();				// Alloca nuova Matrix vuota
				int ir, ic, cc;
				if ((_row > 0) && (m._col > 0))
					{
					DATA sum;
					tmp->dat = new DATA[_row * m._col];				// Alloca spazio.
					if (tmp->dat)
						{
						for (ir = 0; ir < _row; ir++)
							for (ic = 0; ic < m._col; ic++)
								{
								sum = (DATA)0;
								for (cc = 0; cc < _col; cc++)		// Ciclo per sommatoria
									{
									sum += dat[ir * _col + cc] * m.dat[cc * m._col + ic];
									}
								tmp->dat[ir * m._col + ic] = sum;
								}
						tmp->_row = _row;
						tmp->_col = m._col;
						} else
						{
						tmp->_row = tmp->_col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					} else
					{												// Matrice nulla
					tmp->_row = tmp->_col = 0;
					tmp->dat = (DATA*)nullptr;
					}
					#ifdef _DEBUG
					cout << "Matrix::operator*(const Matrix <DATA> &m)" << endl;
					#endif
					return *tmp;
				}
			Matrix <DATA> &operator*=(const Matrix <DATA> &m) requires RQsumdifprod<DATA>
				{
				if (_col != m._row)
					{
					throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
					}
				Matrix<DATA> *tmp = new Matrix<DATA>();				// Alloca nuova Matrix vuota
				int ir, ic, cc;
				if ((_row > 0) && (m._col > 0))
					{
					DATA sum;
					tmp->dat = new DATA[_row * m._col];				// Alloca spazio.
					if (tmp->dat)
						{
						for (ir = 0; ir < _row; ir++)
							for (ic = 0; ic < m._col; ic++)
								{
								sum = (DATA)0;
								for (cc = 0; cc < _col; cc++)		// Ciclo per sommatoria
									{
									sum += dat[ir * _col + cc] * m.dat[cc * m._col + ic];
									}
								tmp->dat[ir * m._col + ic] = sum;
								}
						tmp->_row = _row;
						tmp->_col = m._col;
						} else
						{
						tmp->_row = tmp->_col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					} else
					{												// Matrice nulla
					tmp->_row = tmp->_col = 0;
					tmp->dat = (DATA*)nullptr;
					}
					if (dat != (DATA*)nullptr)							// Dealloca (in this) lo spazio per la matrice 
						delete[] dat;
					dat = tmp->dat;										// Assegna (a this) i valori di tmp 
					_row = tmp->_row;
					_col = tmp->_col;
					#ifdef _DEBUG
					cout << "Matrix::operator*=(const Matrix <DATA> &m)" << endl;
					#endif
					return *this;
				}
			// Prodotto scalare (solo vettori)
			DATA &operator^(const Matrix <DATA> &m) requires RQsumdifprod<DATA>
				{
				if ((_row != m._row) || (_col != m._col) || (_row == 0) || (_col == 0))
					{
					throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
					}
				if ((_row != 1) && (_col != 1))
					{
					throw std::runtime_error(MatrixDef::ERR_NOTVECTOR);
					}
				DATA *tmp = new DATA;
				*tmp = (DATA)0;
				int ir, ic;
				for (ir = 0; ir < _row; ir++)
					for (ic = 0; ic < _col; ic++)
						*tmp += dat[ir * _col + ic] * m.dat[ir * _col + ic];
				#ifdef _DEBUG
				cout << "Matrix::operator^(const Matrix <DATA> &m)" << endl;
				#endif
				return *tmp;
				}
			// Prodotto con scalare (T * Matrix<T> e Matrix<T> * T) e divisione (solo Matrix<T> / T).
			Matrix <DATA> &operator*(const DATA &x) requires RQsumdifprod<DATA>
				{
				Matrix<DATA> *tmp = new Matrix<DATA>();							// Alloca nuova Matrix vuota
				int ir, ic;
				if ((_row > 0) && (_col > 0))
					{
					tmp->dat = new DATA[_row * _col];							// Alloca spazio.
					if (tmp->dat)
						{
						for (ir = 0; ir < _row; ir++)							// Somma
							for (ic = 0; ic < _col; ic++)
								tmp->dat[ir * _col + ic] = dat[ir * _col + ic] * x;
						tmp->_row = _row;
						tmp->_col = _col;
						} else
						{
						tmp->_row = tmp->_col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					} else
					{															// Matrice nulla
					tmp->_row = tmp->_col = 0;
					tmp->dat = (DATA*)nullptr;
					}
					#ifdef _DEBUG
					cout << "Matrix::operator*(const DATA &x)" << endl;
					#endif
					return *tmp;
				}
			friend Matrix<DATA> &operator*(const DATA &sx, const Matrix <DATA> &dx) requires RQsumdifprod<DATA>
				{																// Inline, per evitare errore linker
				Matrix<DATA> *tmp = new Matrix<DATA>();							// Alloca nuova Matrix vuota
				int ir, ic;
				if ((dx._row > 0) && (dx._col > 0))
					{
					tmp->dat = new DATA[dx._row * dx._col];							// Alloca spazio.
					if (tmp->dat)
						{
						for (ir = 0; ir < dx._row; ir++)							// Somma
							for (ic = 0; ic < dx._col; ic++)
								tmp->dat[ir * dx._col + ic] = dx.dat[ir * dx._col + ic] * sx;
						tmp->_row = dx._row;
						tmp->_col = dx._col;
						}
					else
						{
						tmp->_row = tmp->_col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					}
				else
					{															// Matrice nulla
					tmp->_row = tmp->_col = 0;
					tmp->dat = (DATA*)nullptr;
					}
				#ifdef _DEBUG
				cout << "Matrix::operator*(const DATA &sx, const Matrix <DATA> &dx)" << endl;
				#endif
				return *tmp;
				}
			Matrix <DATA> &operator/(const DATA &x) requires RQsumdifprod<DATA>
				{
				Matrix<DATA> *tmp = new Matrix<DATA>();							// Alloca nuova Matrix vuota
				int ir, ic;
				if ((_row > 0) && (_col > 0))
					{
					tmp->dat = new DATA[_row * _col];							// Alloca spazio.
					if (tmp->dat)
						{
						for (ir = 0; ir < _row; ir++)							// Somma
							for (ic = 0; ic < _col; ic++)
								tmp->dat[ir * _col + ic] = dat[ir * _col + ic] / x;		// Possibile throw
						tmp->_row = _row;
						tmp->_col = _col;
						} else
						{
						tmp->_row = tmp->_col = 0;
						throw std::runtime_error(MatrixDef::ERR_ALLOC);
						}
					} else
					{															// Matrice nulla
					tmp->_row = tmp->_col = 0;
					tmp->dat = (DATA*)nullptr;
					}
					#ifdef _DEBUG
					cout << "Matrix::operator/(const DATA &x)" << endl;
					#endif
					return *tmp;
				}
			// Speciali 
			static const Matrix <DATA> &Id(int sz) requires RQsumdifprod<DATA>
				{
				if (sz <= 0)
					{
					throw std::runtime_error(MatrixDef::ERR_WRONGPARAM);
					}
				Matrix<DATA> *tmp = new Matrix<DATA>(sz, sz, (DATA)0);
				int irc;
				for (irc = 0; irc < sz; irc++)
					tmp->dat[irc * sz + irc] = (DATA)1;
				return *tmp;
				}
			void setId(int sz) requires RQsumdifprod<DATA>
			{
				if (sz <= 0)
				{
					throw std::runtime_error(MatrixDef::ERR_WRONGPARAM);
				}
				clear();
				_row = _col = sz;
				dat = new DATA[sz * sz];
				if (!dat)
				{
					_row = _col = 0;
					throw std::runtime_error(MatrixDef::ERR_ALLOC);
				}
				else
				{
					for (int ir = 0; ir < _row; ir++)
						for (int ic = 0; ic < _col; ic++)
						{
							dat[ir * _col + ic] = (ir == ic) ? (DATA)1 : (DATA)0;
						}
				}

#ifdef _DEBUG
			cout << "setId(int sz)" << endl;
#endif
			}
			// Friend Iterator
			#ifdef SPECIALITERATOR
			friend MatrixIterator<DATA>;
			int get_iterators_num()
				{
				return _iterators;
				}
			#endif
			iterator begin() noexcept { return iterator(dat); }
			iterator end() noexcept { return iterator(dat+(_row*_col)); }

			// Friend ostream operators
			template<typename DATA> friend ostream& operator<<(ostream& os, Matrix<DATA>& m);
			template<typename DATA> friend ostream& operator<<(ostream& os, const Matrix<DATA>& m);
		};

	template <typename DATA> DATA *Matrix<DATA>::_empty = new DATA;				// Costruttore di default di DATA. Deallocato a fine programma.
	template <typename DATA> size_t Matrix<DATA>::_datasize = sizeof(*_empty);	// Dimensione di DATA

	template<typename U> ostream& operator<<(ostream& ostr, Matrix<U>& m)
		{
		ostr << m.to_string();
		return ostr;
		}
	template<typename U> ostream& operator<<(ostream& ostr, const Matrix<U>& m)
	{
		ostr << m.to_string();
		return ostr;
	}

	namespace linearsys		// nested
		{
		template <typename DATA, class MOD> class LinearSys
			{
			private:
				Matrix<DATA> a;				// Matrice fattorizzata
				Matrix<int> pivot;			// Vettore di pivot
				//const MOD EPS = std::numeric_limits<MOD>::epsilon();	// Superfluo. Usare std::numeric_limits<> all'interno delle funzioni (con controllo di requires<>).
				//const MOD MINV = std::numeric_limits<MOD>::min();		// Errore, probabilmente rif. ambiguo.

				MOD _epszero;				// Valore minimo considerato > 0
				MOD _det;					// Determinante

			public:

				/* Ctor. Vincolo nel costruttore, non serve altrove */
				LinearSys() requires RQfloatabs<DATA, MOD>
					{
					_epszero = std::numeric_limits<MOD>::epsilon();
					_det = (MOD)0;
					#ifdef _DEBUG		
					cout << "LinearSys(): epsilon=" << _epszero << endl;
					#endif
					}

				/* Soglia sotto la quale un valore si considera nullo */
				void set_eps_zero(MOD eps_zero)
					{
					if (eps_zero > std::numeric_limits<MOD>::epsilon())
						{
						_epszero = eps_zero;
						}
					}
				void set_eps_zero_ratio(MOD eps_zero)
					{
					if (eps_zero >= (MOD)1)
						{
						_epszero = eps_zero * std::numeric_limits<MOD>::epsilon();
						}
					}
				MOD get_eps_zero()
					{
					return _epszero;
					}

				/* to_string() */
				string to_string()
					{
					stringstream ss;
					ss << "Eps_zero=" << _epszero << endl;
					ss << "Det=" << _det << endl;
					ss << "PLU=" << a.to_string() << endl;
					ss << "Pivot=" << pivot.to_string() << endl;
					return ss.str();
					}

				/* Fattorizzazione LU con pivoting parziale e soluzione. Monegato Metodi e algoritmi... CLUT 2008, pag. 41 e succ.*/
				bool factor(Matrix <DATA> &A)
					{
					int n;
					n = A.rows();
					if (n != A.cols())
						{
						throw std::runtime_error(MatrixDef::ERR_NOTSQUARE);
						}
					if (n < 1)
						{
						throw std::runtime_error(MatrixDef::ERR_ZEROSIZE);
						}
					pivot.dim(n - 1, 1);		// Vettore colonna n-1 elementi
					a = A;							// Assegnazione della classe Matrix<>
					if (n == 1)						// Se matrice 1x1: il determinante è il valore.
						{
						_det = a(0, 0);
						if (abs(_det) < _epszero)
							{
							return false;
							}
						return true;
						}
					int k, i, io, j;				// Ciclo di calcolo
					MOD amax;
					DATA tmp;
					for (_det = (DATA)1, k = 0; k < n - 1; k++)			// Ciclo 1 su tutte le righe k=0...n-2
						{
						for (amax = (MOD)0, io = i = k; i < n; i++)			// Cerca la riga (io), dalla k in poi, con il massimo amax sulla colonna k.
							{
							if (abs(a(i, k)) >= amax)
								{
								io = i;
								amax = abs(a(i, k));
								}
							}
						pivot(k, 0) = io;				// Imposta il pivot
						if (amax < _epszero)					// Se il valore massimo dei coefficienti è (quasi) zero...
							{
							_det = (MOD)0;						// ...azzera il determinante ed esce.
							return false;
							}
						if (io != k)							// Se l'indice non è k...
							{
							for (j = k; j < n; j++)				// ...scambia le righe k e io (solo la parte superiore, l'inferiore è nulla)...
								{
								tmp = a(k, j);
								a(k, j) = a(io, j);
								a(io, j) = tmp;
								}
							_det = -_det;							// ...e cambia segno al determinante.
							}
						for (i = k + 1; i < n; i++)				// Ciclo 2, per tutte le righe successive alla riga k attuale
							{
							a(i, k) = -a(i, k) / a(k, k);		// Coefficiente per metodo di Gauss (su riga i, colonna k), sotto diagonale.
							for (j = k + 1; j < n; j++)			// Combina linearmente le due righe, sopra la diagonale 
								{
								a(i, j) = a(i, j) + a(i, k) * a(k, j);
								}
							}								// Fine ciclo 2
						_det = _det * a(k, k);
						}
					if (abs(a(n - 1, n - 1)) < _epszero)
						{
						_det = (MOD)0;
						return false;
						}
					_det = _det * a(n - 1, n - 1);
					return true;
					}
				bool solve_check(bool throw_exception = true)
					{
					bool ok = true;
					int n = a.rows();
					if (n != a.cols())
						{
						ok = false;
						if (throw_exception)		throw std::runtime_error(MatrixDef::ERR_NOTSQUARE);
						}
					if (n < 1)
						{
						ok = false;
						if (throw_exception)	throw std::runtime_error(MatrixDef::ERR_ZEROSIZE);
						}
					if (abs(_det) < _epszero)
						{
						ok = false;
						if (throw_exception)	throw std::runtime_error(MatrixDef::ERR_SINGULAR);
						}
					if ((pivot.rows() != n - 1) || (pivot.cols() != 1))
						{
						ok = false;
						if (throw_exception)	throw std::runtime_error(MatrixDef::ERR_PIVOT);
						}
					return ok;
					}
				bool solve(Matrix <DATA> &x, Matrix <DATA> &b)
					{
					int n = a.rows();
					solve_check();
					if ((b.rows() != n) || (b.cols() != 1))		// Verifica dimensioni vettore termini noti
						{
						throw std::runtime_error(MatrixDef::ERR_SZMISMATCH);
						}
					x = b;
					if (n == 1)
						{
						if (abs(a(0, 0)) < _epszero)
							{
							throw std::runtime_error(MatrixDef::ERR_SINGULAR);
							}
						x(0, 0) = x(0, 0) / a(0, 0);
						return true;
						}
					int k, j, i;
					DATA tmp;

					for (k = 0; k < n - 1; k++)			// Ciclo su tutte le righe, tranne l'ultima
						{
						j = pivot(k, 0);		// Ottiene la riga da scambiare con le riga k
						if (j != k)						// Le scambia, se necessario
							{
							tmp = x(j, 0);
							x(j, 0) = x(k, 0);
							x(k, 0) = tmp;
							}
						for (i = k + 1; i < n; i++)
							{
							x(i, 0) += a(i, k) * x(k, 0);
							}
						}
					x(n - 1, 0) = x(n - 1, 0) / a(n - 1, n - 1);
					for (i = n - 2; i >= 0; i--)
						{
						for (tmp = (DATA)0.0, j = i + 1; j < n; j++)
							tmp += a(i, j) * x(j, 0);
						x(i, 0) = (x(i, 0) - tmp) / a(i, i);
						}
					return true;
					}

			};


		}	// end of namespace linearsys

	}	// end of namespace matrix
#endif