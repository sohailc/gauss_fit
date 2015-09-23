#ifndef GAUSS_NEWTON
#define GAUSS_NEWTON

#include "linearLeastSquares.h"

int gaussNewton(matrix x, 
				matrix y, 
				matrix paramsInit,
				int (*fitFunction)(matrix,matrix,matrix*),
				double (*derivatives[])(matrix, matrix),
				matrix *paramsFinal);
				
#endif
