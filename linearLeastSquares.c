/* This code solves A matrix equation in the form [X][b] = [Y]
the matrix [b] is n-by-o, [X] is a m-by-n matrix and [Y] is a m-by-o matrix
here we demand m >= n 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "linearLeastSquares.h"

void printMat(matrix A)
{
    int i, j;
    double toPrint;
    
    for (i=0;i<A.m;i++)
    {
        for (j=0;j<A.n;j++)
        {
            toPrint = EL(A,i,j);
            printf("%.4le ",toPrint);
        }
        printf("\n");
    }
    
}

void fprintMat(FILE *fh, matrix A)
{
    int i, j;
    double toPrint;
    
    for (i=0;i<A.m;i++)
    {
        for (j=0;j<A.n;j++)
        {
            toPrint = EL(A,i,j);
            fprintf(fh, "%.4le ",toPrint);
        }
        fprintf(fh, "\n");
    }
    
}

double vectorMagnitude(matrix A)
{
	double sqrSum = 0, ai;
	int i;
	
	for(i=0;i<A.m;i++)
	{
		ai = EL(A,i,0);
		sqrSum += pow(ai,2);
	}
	
	double result = sqrt(sqrSum);
	
	return result;
}

int elementWisePower(matrix A, double power, matrix *Apower)
{   
    Apower->m = A.m;
	Apower->n = A.n;
	NEW_MAT((*Apower));
    
    int i;
    
    for (i=0;i<A.n*A.m;i++)
        Apower->data[i] = pow(A.data[i], power);
	
    return EXIT_SUCCESS;
}

int matrixAppendColumns(matrix A, matrix B, matrix *AB)
{
    if (A.m != B.m)
    {
	printf("Error: append matrixes column wise: the number of row nust be the same\n");
	return EXIT_FAIL;
    }

    AB->m = A.m;
    AB->n = A.n + B.n;
    NEW_MAT((*AB));
    
    int i, j;
    
    for (i=0;i<A.m;i++)
    {
	for (j=0;j<A.n+B.n;j++)
	{
	    if (j<A.n)
	    {
		EL((*AB),i,j) = EL(A,i,j);
	    }else{
		EL((*AB),i,j) = EL(B,i,j-A.n);
	    }
	}
    }
    
    return EXIT_SUCCESS;
}

int matrixCopy(matrix A, matrix *B)
{
	B->m = A.m;
	B->n = A.n;
	NEW_MAT((*B));
	
	int i, j;
	
	for (i=0;i<A.m;i++)
	{
		for (j=0;j<A.n;j++)
		{
			EL((*B),i,j) = EL(A,i,j);
		}
	}
	
	return EXIT_SUCCESS;
}

int matrixGetColumn(matrix A, int column, matrix *colMat)
{
    if (column > A.n)
    {
		printf("Error: get a column from a matrix; column to get is larger then the number of columns in the matrix!\n");
		return EXIT_FAIL;
    }
    
    colMat->m = A.m;
    colMat->n=1;
    NEW_MAT((*colMat));
    
    int i;
    
    for (i=0;i<A.m;i++)
    {
		EL((*colMat),i,0) = EL(A,i,column);
    }
    
    return EXIT_SUCCESS;
}

int matrixGetRow(matrix A, int row, matrix *rowMat)
{
    if (row > A.m)
    {
		printf("Error: get a row from a matrix; row index to get is larger then the number of row in the matrix!\n");
		return EXIT_FAIL;
    }
    
    rowMat->m = 1;
    rowMat->n = A.n;
    NEW_MAT((*rowMat));
    
    int i;
    
    for (i=0;i<A.n;i++)
    {
		EL((*rowMat),0,i) = EL(A,row,i);
    }
    
    return EXIT_SUCCESS;
}

int appendPowerColumnsToMat(matrix A, int maxPower, bool beginWithPowerZero, matrix *Apower)
{
    
    // Case: beginWithPowerZero = false
    // if A=[x1,x2] and maxPower = 3, this will return [x1,x2,x1^2,x2^2,x1^3,x2^3]
    // Case: beginWithPowerZero = true
    // if A=[x1,x2] and maxPower = 3, this will return [1,1,x1,x2,x1^2,x2^2,x1^3,x2^3]
    
    if (beginWithPowerZero == false)
    {
    	Apower->n = A.n*maxPower;
    }else{
    	Apower->n = A.n*maxPower+1;
    }
    
    Apower->m = A.m;
    NEW_MAT((*Apower));
    
    int i, j, power;
    
    for (i=0;i<A.m;i++)
    {
		for (j=0;j<A.n;j++) // and NOT (j=0;j<A.n*maxPower;j++)
		{
			if (beginWithPowerZero == false)
			{
				for (power=0;power<maxPower;power++)
				{
					EL((*Apower),i,j+power*A.n) = pow( EL(A,i,j), power+1);
				}
			}else{
				for (power=0;power<=maxPower;power++)
				{
					EL((*Apower),i,j+power*A.n) = pow( EL(A,i,j), power);
				}
			}
		}
    }
    
    
    return EXIT_SUCCESS;
    
}

int matrixDifference(matrix A, matrix B, matrix *D)
{
	if ((A.m != B.m) || A.n != B.n)
	{
		printf("Error: matrix difference: matrix dimensions must be the same");
		return EXIT_FAIL;
	}
	
	D->m = A.m;
	D->n = A.n;
	NEW_MAT((*D));

	double v1, v2;
	int i;
	int j;
	
	for (i=0;i<A.m;i++)
	{
		for (j=0;j<A.n;j++)
		{
			v1 = EL(A,i,j);
			v2 = EL(B,i,j);
			EL((*D),i,j) = v1-v2;
		}
	}
	
	return EXIT_SUCCESS;
}

int matrixMultiplication(matrix A, matrix B, matrix *AtimesB)
{
    if (A.n != B.m)
    {
	printf("Error: matrix multiplication: the number of columns of the left matrix must be the same as the number of row in the right matrix\n");
	return EXIT_FAIL;
    }
    
    AtimesB->m = A.m;
    AtimesB->n = B.n;
    NEW_MAT((*AtimesB));
    
    int i, j, k;
    double sum;
    
    for(i=0;i<A.m;i++) // a row of A
    {
	for (j=0;j<B.n;j++) // a column of B
	{
	    sum = 0;
	    for (k=0;k<A.n;k++)
	    {
		sum += EL(A,i,k) * EL(B,k,j);
	    }
	    EL((*AtimesB),i,j) = sum;
	}
    }
    
    return EXIT_SUCCESS;
}

int choleskyDecomposition(matrix X, matrix *L) // A MUST be a square matrix of size shape-by-shape
                                                            
{
    if (X.m != X.n)
    {
	printf("Error: cholesky decomposition; matrix must be square\n");
	return EXIT_FAIL;
    }

    int i,j,k; // loop variables
    
    L->m = X.m;
    L->n = X.n;
    NEW_MAT((*L));
    double tmp;

    double sum = 0.0;
    
    for (i=0;i<X.n;i++)
    {
        for (k=0;k<i+1;k++)
			sum += pow( EL((*L),i,k), 2);
        
		tmp = EL(X,i,i)-sum;
        
        if (tmp < 0)
        {
            printf("Error: Matrix not positive definite\n");
            return EXIT_FAIL;
        }
        
        EL((*L),i,i) = sqrt(tmp);
        
        sum = 0.0;
        
        for(j=i+1;j<X.n;j++)
        {
            for(k=0;k<i+1;k++)
				sum += EL((*L),i,k) * EL((*L),j,k);
            
            EL((*L),j,i) = (EL(X,j,i)-sum)/EL((*L),i,i);
            
            sum = 0.0;
        }
    }
    
    return EXIT_SUCCESS;
}

int matrixTranspose(matrix A, matrix *AT)
{
    AT->m=A.n;
    AT->n=A.m;
    NEW_MAT((*AT));
    
    int i,j;
    
    for(i=0;i<A.m;i++)
    {
		for(j=0;j<A.n;j++)
		{
			EL((*AT),j,i) = EL(A,i,j);
		}
    }
    
    return EXIT_SUCCESS;
}

int ATransposeTimesB(matrix A, matrix B, matrix *ATB)
{
    if (A.m != B.m)
    {
		printf("Error: A^T * B; the number of rows in A should be the same as the number of rows in B");
		return EXIT_FAIL;
		return EXIT_FAIL;
    }
    
    ATB->m = A.n;
    ATB->n = B.n;
    NEW_MAT((*ATB));

    int i, j, k;
    double sum = 0;
    
    for (i=0;i<A.n;i++)
    {
        for (j=0;j<B.n;j++)
        {
            sum = 0;
            for (k=0;k<A.m;k++)
			sum += EL(A,k,i) * EL(B,k,j); // we multiply column i with column j

			EL((*ATB),i,j) = sum;
            
        }
    }

    return EXIT_SUCCESS;
}

int ATimesBTranspose(matrix A, matrix B, matrix *ABT)
{
    if (A.n != B.n)
    {
		printf("Error: A * B^T; the number of column in A should be the same as the number of columns in B");
		return EXIT_FAIL;
    }
    
    ABT->m = A.m;
    ABT->n = B.m;
    NEW_MAT((*ABT));
    
    int i, j, k;
    double sum = 0;
    
    for (i=0;i<A.m;i++)
    {
        for (j=0;j<B.m;j++)
        {
            sum = 0;
            for (k=0;k<A.n;k++)
				sum += EL(A,i,k) * EL(B,j,k); // we multiply column i with column j

			EL((*ABT),i,j) = sum;
            
        }
    }
    
    return EXIT_SUCCESS;
}

int performBackSubstitution(matrix L, matrix u, matrix *w) // given a upper triangle matrix [L^T], solve [L^T][w] = [u]
                                                           // Since our implementation of the cholesky decomposition
                                                           // only calculates the lower triangle matrix [L], we will
                                                           // perform the calculation here by simply reversing the
                                                           // indexes in [L]. THUS HERE [L] IS A LOWER TRIANGLE MATRIX
{
    if (L.m != u.m)
    {
        printf("Error: back substitution; the number of rows in L and u need to be the same\n");
        return EXIT_FAIL;
    }
    
    w->m = L.n;
    w->n = u.n;
    NEW_MAT((*w));
    
    int i, j, k;
    double sum;
    
    for (k=0;k<u.n;k++)
    {
        for (i=L.m-1;i>=0;i--)
        {
            sum = 0.0;
            for (j=L.n-1;i<j;j--)
            {
                sum += EL(L,j,i) * EL((*w),j,k);
            }
            EL((*w),i,k) = (EL(u,i,k) - sum)/EL(L,i,i);
        }
    }
    
    return EXIT_SUCCESS;
}

int performForwardSubstitution(matrix L, matrix u, matrix *w) // given a lower triangle matrix [L], solve [L][w] = [u] 
{
    if (L.m != u.m)
    {
        printf("Error: forward substitution; the number of rows in L and u need to be the same\n");
        return EXIT_FAIL;
    }
    
    w->m = L.n;
    w->n = u.n;
    NEW_MAT((*w));
    
    int i, j, k;
    double sum;
    
    for (k=0;k<u.n;k++)
    {
        for (i=0;i<L.m;i++)
        {
            sum = 0.0;
            for (j=0;j<i;j++)
            {
                sum += EL(L,i,j) * EL((*w),j,k);
            }
            EL((*w),i,k) = (EL(u,i,k) - sum)/EL(L,i,i);
        }
    }
    
    return EXIT_SUCCESS;
}

int linearLeastSquares(matrix X, matrix Y, matrix *b)
{
    // Solve [X][b] = [Y], the matrix [b] is n-by-o, [X] is a m-by-n matrix and [Y] is a m-by-o matrix
    // The normal equation is [X^T][X][b] = [X^T][Y]. With a cholesky factorization, 
    // this reduces to [L][L^T][b] = [X^T][Y]. This is easier since L is a lower triangle matrix
    // We define [L^T][b] = [w] and solve [L][w] = [X^T][Y] with forward-substitution. We then solve
    // [L^T][b] = [w] for [b] with back-subsitution
    
    if (X.n > X.m)
    {
		printf("Error: the number of unknowns may not be greater then the number of equations\n");
		return EXIT_FAIL;
    }
    
    if (X.m != Y.m)
    {
        printf("Error: the number of data points left and right must be the same\n");
        return EXIT_FAIL;
    }

    // lets first compute [X^T][X] and [X^T][Y]
    matrix XTtimesX, XTtimesY;
    EVAL(ATransposeTimesB(X,X,&XTtimesX));
    EVAL(ATransposeTimesB(X,Y,&XTtimesY));

    // lets then calculate the cholesky decomposition of [X^T][X]
    matrix L;
    EVAL(choleskyDecomposition(XTtimesX,&L));
    
    matrix w;
    EVAL(performForwardSubstitution(L,XTtimesY,&w));
    EVAL(performBackSubstitution(L,w,b));
    
    free(XTtimesX.data);
    free(XTtimesY.data);
    free(L.data);
    free(w.data);
    
    return EXIT_SUCCESS;
    
}

int polyFit(matrix x, matrix y, int order, bool throughZero, matrix *params)
{
    if (x.n!=1 || y.n!=1)
    {
	printf("Error: polyfit; one dimensional data expected\n");
	return EXIT_FAIL;
    }
    
    bool appendZeroOrder;
    if (throughZero)
        appendZeroOrder = false;
    else
        appendZeroOrder = true;
    
    matrix xPowers;
    EVAL(appendPowerColumnsToMat(x,order,appendZeroOrder,&xPowers));
    EVAL(linearLeastSquares(xPowers,y,params));
    
    free(xPowers.data);
    
    return EXIT_SUCCESS;
    
}

