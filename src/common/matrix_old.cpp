/*
 * matrix.cpp
 *
 *  Created on: 22 nov. 2013
 *      Author: jfellus
 */

#include "matrix.h"


static Matrix* cur_cov_mt_C;
static Matrix* cur_cov_mt_X;
static __multithread__(compute_cov_mt) (int i) {
	Matrix& X = *cur_cov_mt_X;
	Matrix& C = *cur_cov_mt_C;

	for(int j=i; j<X.width; j++) {
		C[i+j*X.width] = 0;
		for(int k=0; k<X.height; k++) C[i+j*X.width] += X[i + k*X.width] * X[j + k*X.width];
		C[i+j*X.width] /= X.height;

		C[j+i*X.width] = C[i+j*X.width];
	}
}


void Matrix::covariance_mt(Matrix& X) {
	cur_cov_mt_C = this;
	cur_cov_mt_X = &X;
	compute_cov_mt(X.width);
}

