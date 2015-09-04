/*
 * math.cpp
 *
 *  Created on: 21 janv. 2014
 *      Author: jfellus
 */

#include "math.h"
#include "matrix.h"
#include "matrix_double.h"



float frand() {
	return ((float)rand())/INT_MAX;
}

float frand(float min, float max) {
	return frand()*(max-min)+min;
}

void randvec(float* v, int n, float min, float max) {
	for(int i=0; i<n; i++) {
		v[i] = frand(min,max);
	}
}

void randvec(double* v, int n, double min, double max) {
	for(int i=0; i<n; i++) {
		v[i] = frand(min,max);
	}
}


float rand_gaussian(float mu, float sigma) {
	float x1, x2, w, y1;

	do {
		x1 = 2.0 * frand() - 1.0;
		x2 = 2.0 * frand() - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;
	return (y1*sigma) + mu;
}


float rand_exp(float lambda) {
	float u = frand();
	return -log(u)/lambda;
}




void randvec_gaussian(float* v, int d) {
	for(int i=0; i<d; i++) v[i] = rand_gaussian(0,1);
}

static MatrixFloat _randvec_gaussian_tmp;
void randvec_gaussian(float* v, MatrixFloat& mu, MatrixFloat& sigma) {
//	if(_randvec_gaussian_tmp.width != mu.width)
//		_randvec_gaussian_tmp.realloc(mu.width, 1);

	randvec_gaussian(_randvec_gaussian_tmp, mu.width);
	for(int i=0; i<mu.width; i++) v[i] = mu[i];
	matrix_CpAtB_float(v, sigma, _randvec_gaussian_tmp, sigma.width, sigma.width, 1);
}

void rand_covariance(MatrixFloat& cov, float s, int rank) {
	Matrix x(cov.width, rank);
	x.rand(-s,s);
	cov.covariance(x);
}


void randvec_gaussian(double* v, int d) {
	for(int i=0; i<d; i++) v[i] = rand_gaussian(0,1);
}

static MatrixDouble _randvec_gaussian_tmp_d;
void randvec_gaussian(double* v, MatrixDouble& mu, MatrixDouble& sigma) {
	if(_randvec_gaussian_tmp_d.width != mu.width)
		_randvec_gaussian_tmp_d.realloc(mu.width, 1);

	randvec_gaussian(_randvec_gaussian_tmp_d, mu.width);
	for(int i=0; i<mu.width; i++) v[i] = mu[i];
	matrix_CpAtB_double(v, sigma, _randvec_gaussian_tmp_d, sigma.width, sigma.width, 1);
}

void rand_covariance(MatrixDouble& cov, double s, int rank) {
	MatrixDouble x(cov.width, rank);
	x.rand(-s,s);
	cov.covariance(x);
}
