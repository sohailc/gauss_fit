#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "gaussNewton.h"
#include "linearLeastSquares.h"

int fitFunction(matrix xy, matrix params, matrix *modelValues)
{
	/*
	
	The function z = A*exp(-(x-mux)^2/(2*sigx^2) -  (y-muy)^2/(2*sigy^2)) + z0
	We will use z = A*exp(-wx^2 - wy^2) + z0, with wx = (x-mux)/(sqrt(2)*sigx), wy = (y-muy)/(sqrt(2)*sigy)
	
	*/

	int M = xy.m;
	int i;
	
	double mux = EL(params,0,0);
	double muy = EL(params,1,0);
	
	double sigx = EL(params,2,0);
	double sigy = EL(params,3,0);
	
	double A = EL(params,4,0);
	double z0 = EL(params,5,0);
	double xi, yi, zi;
	
	double wx, wy;
	
	for(i=0;i<M;i++)
	{
		xi = EL(xy,i,0);
		yi = EL(xy,i,1);
		
		wx = (xi-mux)/(sqrt(2)*sigx);
		wy = (yi-muy)/(sqrt(2)*sigy);
		
		zi = A*exp(-pow(wx,2) - pow(wy,2)) + z0;
		
		EL((*modelValues),i,0) = zi;
	}
	
	return EXIT_SUCCESS;
}

double derivative_ROverMux(matrix xy, matrix params)
{
	// differentiate 
	// z = A*exp(-wx^2 - wy^2) + z0, with wx = (x-mux)/(sqrt(2)*sigx), wy = (y-muy)/(sqrt(2)*sigy)
	// with respect to mux (we return the negative derivative -dz/dmux !)
	// z = A*exp(-wx^2)*exp(-wy^2) + z0
	// dz/dmux = dz/dwx * dwx/dmux = A*exp(-wx^2)*exp(-wy^2)*-2*wx * -1/(sqrt(2)*sigx) = 2*wx*A*exp(-wx^2)*exp(-wy^2)/(sqrt(2)*sigx)
	
	double xi = EL(xy,0,0), yi = EL(xy,0,1);
	
	double mux = EL(params,0,0);
	double muy = EL(params,1,0);
	
	double sigx = EL(params,2,0);
	double sigy = EL(params,3,0);
	
	double A = EL(params,4,0);
	
	double wx = (xi-mux)/(sqrt(2)*sigx);
	double wy = (yi-muy)/(sqrt(2)*sigy);
	
	return -2*wx*A*exp(-pow(wx,2)-pow(wy,2))/(sqrt(2)*sigx);
}

double derivative_ROverMuy(matrix xy, matrix params)
{
	// differentiate 
	// z = A*exp(-wx^2 - wy^2) + z0, with wx = (x-mux)/(sqrt(2)*sigx), wy = (y-muy)/(sqrt(2)*sigy)
	// with respect to muy (we return the negative derivative -dz/dmuy !)
	// z = A*exp(-wx^2)*exp(-wy^2) + z0
	// dz/dmuy = dz/dwy * dwy/dmux = A*exp(-wx^2)*exp(-wy^2)*-2*wy * -muy/(sqrt(2)*sigy) = 2*wy*muy*A*exp(-wx^2)*exp(-wy^2)/(sqrt(2)*sigy)
	
	double xi = EL(xy,0,0), yi = EL(xy,0,1);
	
	double mux = EL(params,0,0);
	double muy = EL(params,1,0);
	
	double sigx = EL(params,2,0);
	double sigy = EL(params,3,0);
	
	double A = EL(params,4,0);
	
	double wx = (xi-mux)/(sqrt(2)*sigx);
	double wy = (yi-muy)/(sqrt(2)*sigy);
	
	return -2*wy*A*exp(-pow(wx,2)-pow(wy,2))/(sqrt(2)*sigy);
}

double derivative_ROverSigx(matrix xy, matrix params)
{
	// differentiate 
	// z = A*exp(-wx^2 - wy^2) + z0, with wx = (x-mux)/(sqrt(2)*sigx), wy = (y-muy)/(sqrt(2)*sigy)
	// with respect to sigx (we return the negative derivative -dz/dsigx !)
	// z = A*exp(-wx^2)*exp(-wy^2) + z0
	// dz/dsigx = dz/dwx * dwx/dsigx = A*exp(-wx^2)*exp(-wy^2)*-2*wx * -(x-mux)/(sqrt(2)*sigx^2)

	double xi = EL(xy,0,0), yi = EL(xy,0,1);
	
	double mux = EL(params,0,0);
	double muy = EL(params,1,0);
	
	double sigx = EL(params,2,0);
	double sigy = EL(params,3,0);
	
	double A = EL(params,4,0);
	
	double wx = (xi-mux)/(sqrt(2)*sigx);
	double wy = (yi-muy)/(sqrt(2)*sigy);
	
	return -2*A*wx*(xi-mux)/(sqrt(2)*pow(sigx,2))*exp(-pow(wx,2)-pow(wy,2));
}

double derivative_ROverSigy(matrix xy, matrix params)
{
	// differentiate 
	// z = A*exp(-wx^2 - wy^2) + z0, with wx = (x-mux)/(sqrt(2)*sigx), wy = (y-muy)/(sqrt(2)*sigy)
	// with respect to sigy (we return the negative derivative -dz/dsigy !)
	// z = A*exp(-wx^2)*exp(-wy^2) + z0
	// dz/dsigy = dz/dwy * dwy/dsigy = A*exp(-wx^2)*exp(-wy^2)*-2*wy * -(y-muy)/(sqrt(2)*sigy^2)

	double xi = EL(xy,0,0), yi = EL(xy,0,1);
	
	double mux = EL(params,0,0);
	double muy = EL(params,1,0);
	
	double sigx = EL(params,2,0);
	double sigy = EL(params,3,0);
	
	double A = EL(params,4,0);
	
	double wx = (xi-mux)/(sqrt(2)*sigx);
	double wy = (yi-muy)/(sqrt(2)*sigy);
	
	return -2*A*wy*(yi-muy)/(sqrt(2)*pow(sigy,2))*exp(-pow(wx,2)-pow(wy,2));
}

double derivative_ROverA(matrix xy, matrix params)
{
	// differentiate 
	// z = A*exp(-wx^2 - wy^2) + z0, with wx = (x-mux)/(sqrt(2)*sigx), wy = (y-muy)/(sqrt(2)*sigy)
	// with respect to A (we return the negative derivative -dz/dA !)
	
	double xi = EL(xy,0,0), yi = EL(xy,0,1);
	
	double mux = EL(params,0,0);
	double muy = EL(params,1,0);
	
	double sigx = EL(params,2,0);
	double sigy = EL(params,3,0);
	
	double wx = (xi-mux)/(sqrt(2)*sigx);
	double wy = (yi-muy)/(sqrt(2)*sigy);
	
	return -exp(-pow(wx,2) - pow(wy,2));
}

double derivative_ROverZnull(matrix x, matrix params)
{
	return -1;
}

int estimateParams(matrix xy, matrix z, matrix *estimates)
{
	/*
	An improvement of this function would be to sort the z-values and take the mean of the bottom x% 
	as minZ
	*/

	int i, M = xy.m;
	double mux = 0, muy=0;
	double sigx = 0, sigy = 0;
	double xi, yi, zi, totalZ = 0;
	double maxY = -1E6, minY = 1E6;
	double maxX = -1E6, minX = 1E6;
	double maxZ = -1E6, minZ = 1E6;
	
	for(i=0;i<M;i++)
	{	
		xi = EL(xy,i,0);
		yi = EL(xy,i,1);
		zi = EL(z,i,0);
	
		if (xi > maxX) {maxX = xi;}
		if (xi < minX) {minX = xi;}	
		if (yi > maxY) {maxY = yi;}
		if (yi < minY) {minY = yi;}
		if (zi > maxZ) {maxZ = zi;}
		if (zi < minZ) {minZ = zi;}	
	}
	
	double A = maxZ-minZ;
	
	for(i=0;i<M;i++)
	{
		xi = EL(xy,i,0);
		yi = EL(xy,i,1);
		zi = EL(z,i,0) - minZ;
		
		if (zi < 0.2*A)
			continue; 
			
		totalZ += zi; 
		mux += xi*zi;
		muy += yi*zi;
	}
	
	mux /= totalZ;
	muy /= totalZ;

	double dx = (maxX - minX)/sqrt(M); // we assume that the sampling distance in x and y is approx.
	double dy = (maxY - minY)/sqrt(M); // the same, so the number of samples along x and y is ~sqrt(M)
	
	double volumeUnderGauss = totalZ*dx*dy;
	
	//                       inf inf
	//						 /    /
	// 1/(2*pi*sigx*sigy)  * |    | exp(-wx^2-wy^2) dx dy = 1
	//						 /    /
	//                    -inf  -inf
	
	// A = volumeUnderGauss / (2*pi*sigx*sigy)
	// We will assume sigx = sigy as an initial estimate
	// 2*pi*sigx^2 = volumeUnderGauss / A
	// sig^2 = volumeUnderGauss/A * 1/(2*pi)
	
	sigx = sqrt(volumeUnderGauss/A * 1.0/(2*3.14)); 
	sigy = sigx;
	
	estimates->m = 6;
	estimates->n = 1;
	NEW_MAT((*estimates));
	estimates->data[0] = mux;
	estimates->data[1] = muy;
	estimates->data[2] = sigx;
	estimates->data[3] = sigy;
	estimates->data[4] = A;
	estimates->data[5] = minZ;
	
	return EXIT_SUCCESS;
}

int gaussian2DFit(	double x[], 
					double y[], 
					double intensities[], 
					int dataSize, 
					double (*paramsArray)[])
{
	int result = EXIT_SUCCESS;
	
	matrix xy, z;
	xy.m = dataSize;
	xy.n = 2;
	NEW_MAT(xy);
	
	z.m = dataSize;
	z.n = 1;
	NEW_MAT(z);
	
	int i;
	for(i=0;i<dataSize;i++)
	{
		EL(xy,i,0) = x[i];
		EL(xy,i,1) = y[i];
		EL(z,i,0) = intensities[i];
	}

	matrix paramsInit, params;
	EVAL(estimateParams(xy,z,&paramsInit));
	
	double (*derivatives[])(matrix, matrix) = {	derivative_ROverMux, 
												derivative_ROverMuy, 
												derivative_ROverSigx, 
												derivative_ROverSigy, 
												derivative_ROverA, 
												derivative_ROverZnull};
	
	if (gaussNewton(xy, z, paramsInit, fitFunction, derivatives, &params) == EXIT_FAIL)
	{
		printf("data very noisy or contrast too low (or both)\n");
		result = EXIT_FAIL;
	}else{
	
		for(i=0;i<6;i++){(*paramsArray)[i] = params.data[i];}
	
	}
	
	free(params.data);
	free(paramsInit.data);
	free(xy.data);
	free(z.data);

	return result;

}
