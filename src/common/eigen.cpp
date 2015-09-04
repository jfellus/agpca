/*
 * matrix.cpp
 *
 *  Created on: 22 nov. 2013
 *      Author: jfellus
 */


#include "../../eigen/Eigen/Dense"
#include "../../eigen/Eigen/QR"
using Eigen::MatrixXf;
using Eigen::RowVectorXf;
using Eigen::VectorXf;
using Eigen::Map;
using Eigen::HouseholderQR;
using Eigen::SelfAdjointEigenSolver;
using Eigen::Upper;
#include <iostream>

#include "the_eigen.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// OK MATLAB
void QR(float* A, float* R, size_t n, size_t m) {
	Map<MatrixXf> THIS(A, n, m);
	Map<MatrixXf> RR(R, n, n);
	HouseholderQR<MatrixXf> qr(THIS.transpose());
	MatrixXf thinQ(MatrixXf::Identity(m,n));
	THIS = (qr.householderQ() * thinQ).transpose();
	RR = qr.matrixQR().triangularView<Upper>().transpose();
}

// OK MATLAB
void QR(float* A, float* Q, float* R, size_t n, size_t m) {
	Map<MatrixXf> THIS(A, n,m);
	Map<MatrixXf> QQ(Q, n,m);
	Map<MatrixXf> RR(R, n,n);
	HouseholderQR<MatrixXf> qr(THIS.transpose());
	MatrixXf thinQ(MatrixXf::Identity(m,n));
	QQ = (qr.householderQ() * thinQ).transpose();
	RR = (qr.matrixQR().triangularView<Upper>()).transpose();
}

// OK MATLAB
void qr(float* A, float* Q, float* r, size_t n, size_t m) {
	Map<MatrixXf> THIS(A, n, m);
	Map<MatrixXf> QQ(Q, n, m);
	Map<RowVectorXf> RR(r, n);
	HouseholderQR<MatrixXf> qr(THIS.transpose());
	MatrixXf thinQ(MatrixXf::Identity(m,n));
	QQ = (qr.householderQ() * thinQ).transpose();
	RR = qr.matrixQR().diagonal();
}

// OK MATLAB
void qr(float* A, float * r, size_t n, size_t m) {
	Map<MatrixXf> THIS(A, n, m);
	Map<RowVectorXf> RR(r, n);
	HouseholderQR<MatrixXf> qr(THIS.transpose());
	MatrixXf thinQ(MatrixXf::Identity(m,n));
	THIS = (qr.householderQ() * thinQ).transpose();
	RR = qr.matrixQR().diagonal();
}



struct struct_sort_float
{
    size_t i;
    float val;
};

static int compare_sort_float(const void * a, const void * b){
    float v = ((struct struct_sort_float *)b)->val - ((struct struct_sort_float *)a)->val;
    if (v > 0)
        return 1;
    if (v < 0)
        return -1;
    return 0;
}

static void matrix_sortEig_float_t (float* A,float* d,size_t n) {
	size_t i,j;
	struct struct_sort_float * sSort = (struct struct_sort_float *)malloc(n*sizeof(struct struct_sort_float));
	float * Atemp = (float *)malloc(n*n*sizeof(float));

	memcpy(Atemp, A, n*n*sizeof(float));

	for(i = 0 ; i < n ; i++)
	{
		sSort[i].i = i;
		sSort[i].val = d[i];
	}

	qsort(sSort, n, sizeof(struct struct_sort_float), compare_sort_float);

	// Column vectors !!!!!
	for(i = 0 ; i < n ; i++)
	{
		d[i] = sSort[i].val;
		for(j = 0; j<n; j++) A[n*j + i] = Atemp[n*j + sSort[i].i];
	}

	free(Atemp);
	free(sSort);
}



void eigendecompose(float* A, float* U, float* L, size_t n) {
	Map<MatrixXf> X(A, n,n);
	Map<MatrixXf> UU(U, n,n);
	Map<VectorXf> LL(L, n);
	SelfAdjointEigenSolver<MatrixXf> eig(X);
	LL = eig.eigenvalues();
	UU = eig.eigenvectors().transpose();
	matrix_sortEig_float_t(U,L, n);
}
