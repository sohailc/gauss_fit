all:
	gcc -c -Wall -Werror -g -fPIC linearLeastSquares.c gaussNewton.c gauss2DFit.c 
	gcc -shared -o libGaussFit.so linearLeastSquares.o gaussNewton.o gauss2DFit.o
	rm *.o
