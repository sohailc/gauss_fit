#ifndef LINEAR_LEAST_SQUARES
#define LINEAR_LEAST_SQUARES

#define EXIT_SUCCESS 0
#define EXIT_FAIL 1 

#define bool int
#define true 1
#define false 0

// In this code, matrix data is stored in double arrays in a column wise fashion
// So if array A holds data of a 3-by-4 matrix, this should be interpreted as
//
//      A[0] A[1] A[2]  A[3]
//  A = A[4] A[5] A[6]  A[7]
//      A[8] A[9] A[10] A[11]
//
// Then, A.m = 3 and A.n = 4
typedef struct matrix
{
    double *data;
    int m;
    int n;
}matrix;

#define EL(A,i,j) A.data[i*A.n+j] // a way to access element i,j of a matrix. In the above example, EL(A,1,0) will access A[4]
// A macro to allocate memory for a matrix and checking if the allocation was succesful
#define NEW_MAT(A) A.data = (double*)calloc(A.m*A.n,sizeof(double)); if (A.data==0) {printf("MEMORY ERROR"); return EXIT_FAIL;}
// A macro to see if a function return a failure. Handy to implement a cascading failure
#define EVAL(F) if (F==EXIT_FAIL) {return EXIT_FAIL;}

void printMat(matrix A);
void fprintMat(FILE *fh, matrix A);
int elementWisePower(matrix A, double power, matrix *Apower);
int matrixAppendColumns(matrix A, matrix B, matrix *AB);
int matrixCopy(matrix A, matrix *B);
int matrixGetColumn(matrix A, int column, matrix *colMat);
int matrixGetRow(matrix A, int row, matrix *rowMat);
int appendPowerColumnsToMat(matrix A, int maxPower, bool beginWithPowerZero, matrix *Apower);
int matrixDifference(matrix A, matrix B, matrix *D);
int matrixMultiplication(matrix A, matrix B, matrix *AtimesB);
int choleskyDecomposition(matrix X, matrix *L);
int matrixTranspose(matrix A, matrix *AT);
int ATransposeTimesB(matrix A, matrix B, matrix *ATB);
int ATimesBTranspose(matrix A, matrix B, matrix *ABT);
int performBackSubstitution(matrix L, matrix u, matrix *w);
int performForwardSubstitution(matrix L, matrix u, matrix *w);
int linearLeastSquares(matrix X, matrix Y, matrix *b);
double vectorMagnitude(matrix A);
int polyFit(matrix x, matrix y, int order, bool throughZero, matrix *params);

#endif
