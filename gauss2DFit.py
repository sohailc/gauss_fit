import ctypes
import numpy

exit_success = 0
exit_fail = 1

def main(X,Y,Z):
	"""
	This is an python interface function to the gauss fit library. Given a 
	set of X and Y values which represent points where a gaussian function 
	was evaluated, and the function value Z at those locations, return the 
	fit parameters mux, muy, sigx, sigy, A, flr. We fit the function 

	f(X, Y) = A*exp(-(X-mux)^2/(2*sigx**2) - (Y-muy)^2/(2*sigy**2)) + flr

	Input1: X
	Input2: Y
	Input3: Z
	Output: The fit paramaters as a dictionary

	"""
	
	X = numpy.ravel(X)
	Y = numpy.ravel(Y)
	Z = numpy.ravel(Z)
	N = len(X)
	
	Xc = (N*ctypes.c_double)(*X)
	Yc = (N*ctypes.c_double)(*Y)
	Zc = (N*ctypes.c_double)(*Z)
	paramsc = (6*ctypes.c_double)(*(6*[0.0]))
	
	lib = ctypes.CDLL("./libGaussFit.so")
	fitFunction = lib.gaussian2DFit
	fitFunction.restype = ctypes.c_int
	
	result = fitFunction(Xc, Yc, Zc, N, ctypes.byref(paramsc))
	
	if result == exit_success:
	
		params = numpy.array(paramsc)
		names = ["mux", "muy", "sigx", "sigy", "A", "flr"]
		
		pdict = dict(zip(names, params))
	
	else:
		
		pdict = None
	
	return pdict
