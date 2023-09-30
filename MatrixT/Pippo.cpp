

#include "Pippo.h"


ostream& operator<<(ostream& os, const Pippo &pp)
	{
	os << '{' << pp._p << '}';
	return os;
	}



