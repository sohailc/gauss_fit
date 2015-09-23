#ifndef GAUSS_2D
#define GAUSS_2D

#include "linearLeastSquares.h"

int gaussian2DFit(double x[], double y[], double intensities[], int dataSize, double (*paramsArray)[]);
int fitFunction(matrix xy, matrix params, matrix *modelValues); 
				
#endif
