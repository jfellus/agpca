/*
 * matrix.h
 *
 *  Created on: 15 nov. 2013
 *      Author: jfellus
 */

#ifndef GOSSIP_MATRIX_DOUBLE_H_
#define GOSSIP_MATRIX_DOUBLE_H_

#include "../retin/toolbox/core/SharedMatrixDouble.h"
#include "../retin/toolbox/algebra/matrix_double.h"
#include "../retin/toolbox/algebra/vector_double.h"
#include "common/utils.h"

class MatrixDouble : public shared_matrices::MatrixDouble {
public:

	MatrixDouble():shared_matrices::MatrixDouble(){}
	MatrixDouble(size_t w, size_t h):shared_matrices::MatrixDouble(w,h){this->clear();}

	/////////////
	// METHODS //
	/////////////

	void mean_row(MatrixDouble& mean) {
		matrix_sumt_double(mean, data, width, height);
		vector_sdiv_double(mean, height, width);
	}

	void centering(MatrixDouble& mu_out, MatrixDouble& out) {
		out = *this;
		mean_row(mu_out);
		for(int i=0; i<height; i++) vector_sub_double(out.get_row(i), mu_out, width);
	}

	void covariance(MatrixDouble& X) {
		if(X.width > 0) {covariance_mt(X); return;}
		this->clear();
		matrix_CpAAt_double(data, X, X.width, X.height);
		(*this) /= X.height;
	}

	void covariance_mt(MatrixDouble& X);


	void eigendecompose(MatrixDouble& U, MatrixDouble& L) {
		DBG_START("eigendecomposition");

		_check_("matrix");

		U = (*this); L.clear();
		matrix_eigSym_double(U, L, width);
		matrix_sortEig_double(U, L, width);

		for(int i=0; i<width; i++) if(isnan((float)L.data[i])) {
			DBG("ERROR : eigendecompose produced NAN -> recomputing with a little added delta");
			for(int i=0; i<width; i++) data[i*width+i] += 0.00001f;
			eigendecompose(U,L);
			return;
		}

		DBG_END();
	}

	void reconstruct(MatrixDouble& U, MatrixDouble& L) {
		DBG_START("reconstruction");
		this->clear();
		matrix_CpAdAt_double(data, U, L,width,height);
		DBG_END();
	}

	void _check_(const char* what) {
		for(int i=0; i<width; i++) for(int j=0; j<height; j++) {
			if(isnan((float)data[i*width+j])) FATAL("NAN FAULT : in " << what << " [" << i << "," << j << "]");
		}
	}

	inline double min() { return vector_min_double(data, width*height); }
	inline double max() { return vector_max_double(data, width*height); }
	inline double sum() { return vector_sum_double(data, width*height); }

	inline double n2() { return matrix_n2_double(data,width,height); }
	inline double n1() { return vector_n1_double(data,width*height); }
	inline double l2(MatrixDouble& B) { return matrix_l2_double(data,B,height,width); }

	double l2_normalized(MatrixDouble& B) {
		double norm = n2(); if(norm==0) norm=1;
		return l2(B) / norm;
	}

	int count_nonzeros() {
		int n = 0;
		for(int i=0; i<width*height; i++) if(data[i]!=0) n++;
		return n;
	}

	void normalize_l1() { (*this)/=n1(); }


	void AAt(MatrixDouble& A) {	matrix_CpAAt_double(data,A,A.height,A.width); }


	void set_part(int x, int y, MatrixDouble& mat) {
		for(int l=0; l<mat.height; l++)
			memcpy(get_row(l+y)+x, mat.get_row(l), mat.width*sizeof(double));
	}

	void clear_row(int y) {
		memset(get_row(y), 0, sizeof(double)*width);
	}

	void set(int x, int y, double val) {data[x+y*width]=val;}

	void sadd(double s, MatrixDouble& m) {
		for(int x=0; x<width*height; x++) data[x] += s*m[x];
	}

	int nearest_neighbor(double* v) {
		double d,mind=FLT_MAX;
		int nn = -1;
		for(int k=0; k<height; k++) {
			d = vector_l2p2_double(get_row(k),v,width);
			if(d <= mind) {
				mind = d;
				nn = k;
			}
		}
		return nn;
	}

	double nearest_neighbor_distance(double* v) {
		double d,mind=FLT_MAX;
		for(int k=0; k<height; k++) {
			d = vector_l2p2_double(&data[k*width],v,width);
			if(d <= mind) mind = d;
		}
		return mind;
	}

	void row_sdiv(int row, double s) {
		vector_sdiv_double(get_row(row), s, width);
	}

	void row_smul(int row, double s) {
		vector_smul_double(get_row(row), s, width);
	}

	void row_add(int row, double* v) {
		vector_add_double(get_row(row), v, width);
	}

	void CpAtB(MatrixDouble& A, MatrixDouble& B) {
		matrix_CpAtB_double(data, A, B, width, A.height, height);
	}

	void rand(double min = 0, double max = 1);

	void integrate() {
		for(int i=1; i<width*height; i++)
			data[i] = data[i-1] + data[i];
	}

	///////////////
	// OPERATORS //
	///////////////

	void operator=(const int* i) { *((shared_matrices::MatrixDouble*)this) = i; }
	double operator=(double x) { *((shared_matrices::MatrixDouble*)this) = x; return x;}

	void operator/=(double f) {vector_sdiv_double(data, f, width*height);}
	void operator*=(double f) {vector_smul_double(data, f, width*height);}

	MatrixDouble& operator+=(MatrixDouble & m) { vector_add_double(data, m, width*height); return (*this); }
	MatrixDouble& operator-=(MatrixDouble & m) {
		if(m.height==1) {
			for(int i=0; i<height; i++) vector_sub_double(get_row(i), m, width);
		}
		else vector_sub_double(data, m, width*height);
		return (*this);
	}


	/////////
	// DBG //
	/////////

	void dbg_range() { DBG(min() << " -> " << max()); }
};


#ifdef USE_MATRIX_DOUBLE
	typedef MatrixDouble Matrix;
#endif


#endif /* GOSSIP_MATRIX_DOUBLE_H_ */
