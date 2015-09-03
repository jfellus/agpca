/*
 * matrix.cpp
 *
 *  Created on: 22 nov. 2013
 *      Author: jfellus
 */

#include "matrix_double.h"
#include "math.h"


static MatrixDouble* cur_cov_mt_C;
static MatrixDouble* cur_cov_mt_X;
static __multithread__(compute_cov_mt_d) (int i) {
	MatrixDouble& X = *cur_cov_mt_X;
	MatrixDouble& C = *cur_cov_mt_C;

	for(int j=i; j<X.width; j++) {
		C[i+j*X.width] = 0;
		for(int k=0; k<X.height; k++) C[i+j*X.width] += X[i + k*X.width] * X[j + k*X.width];
		C[i+j*X.width] /= X.height;

		C[j+i*X.width] = C[i+j*X.width];
	}
}


void MatrixDouble::covariance_mt(MatrixDouble& X) {
	cur_cov_mt_C = this;
	cur_cov_mt_X = &X;
	compute_cov_mt_d(X.width);
}

void MatrixDouble::rand(double min, double max) {
	for(int i=0; i<width*height; i++)
		data[i] = frand(min,max);
}
