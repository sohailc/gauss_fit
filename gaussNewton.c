#include <stdlib.h>
#include <stdio.h>
#include "linearLeastSquares.h"
#include <math.h>

int computeJacobian(matrix x, matrix params, double (*derivatives[])(matrix, matrix), matrix *jacobian) 
{
	int M = x.m;
	
	if (jacobian->m != M)
	{
		printf("Error: the first dimension of the jacobian needs to be the same size as the number of sample points");
		return EXIT_FAIL;
	}
	
	if (jacobian->n != params.m)
	{
		printf("Error: the second dimension of the jacobian needs to be the same size as the number of parameters");
		return EXIT_FAIL;
	}
	
	int i, j; // loop veriables
	matrix xi;
	double d;
	
	for (i=0;i<M;i++)
	{
		EVAL(matrixGetRow(x,i, &xi));
		
		for (j=0;j<params.m;j++)
		{
			d = (derivatives[j])(xi, params);
			EL((*jacobian),i,j) = d;
		}
		free(xi.data);
	}

	return EXIT_SUCCESS;
}

int iteration(	matrix x,
						matrix y, 
						matrix *params,
						int (*fitFunction)(matrix,matrix,matrix*),
						double (*derivatives[])(matrix, matrix),
						double *precesion)
{
	
	int N = x.m;
	int nParams = params->m;
	matrix fit;
	fit.m = N;
	fit.n = 1;
	NEW_MAT(fit);
	
	fitFunction(x, *params, &fit);
	matrix r;
	EVAL(matrixDifference(y, fit, &r));
	
		
	matrix jacobian; 
	jacobian.m = N;
	jacobian.n = nParams;
	NEW_MAT(jacobian);
	EVAL(computeJacobian(x, *params, derivatives, &jacobian))

	// now solve jacobian * u = r
	matrix u;
	EVAL(linearLeastSquares(jacobian, r, &u));
	
	*precesion = vectorMagnitude(u);
	
	int i;
	for (i=0;i<params->m;i++){params->data[i] -= u.data[i];}
	
	// free all the meory used
	free(fit.data);
	free(r.data);
	free(jacobian.data);
	free(u.data);
	
	return EXIT_SUCCESS;
}

int gaussNewton(matrix x, 
				matrix y, 
				matrix paramsInit,
				int (*fitFunction)(matrix,matrix,matrix*),
				double (*derivatives[])(matrix, matrix),
				matrix *paramsFinal)
{
	double precision = 0.001;
	int countIter = 0, maxIter = 100;
	double d = 1.0;
	EVAL(matrixCopy(paramsInit, paramsFinal));
	
	int exitCode = EXIT_SUCCESS;
	
	while (d > precision && exitCode == EXIT_SUCCESS)
	{
		countIter += 1;
		if (iteration(x, y, paramsFinal, fitFunction, derivatives, &d) == EXIT_FAIL)
		{
			exitCode = EXIT_FAIL;
			printf("fit failed!\n");
		}
		
		if (countIter == maxIter)
		{
			exitCode = EXIT_FAIL;
			printf("Warning: max iter reached!\n");
		}
	}

	return exitCode;
}
