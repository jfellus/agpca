/*
 * matrix.h
 *
 *  Created on: 15 nov. 2013
 *      Author: jfellus
 */

#ifndef GOSSIP_MATRIX_H_
#define GOSSIP_MATRIX_H_

#include "../retin/toolbox/core/SharedMatrix.h"
#include "../retin/toolbox/algebra/matrix_float.h"
#include "../retin/toolbox/algebra/matrix_double.h"
#include "../retin/toolbox/algebra/vector_float.h"
#include "common/utils.h"


class Matrix : public shared_matrices::Matrix {
public:

	Matrix():shared_matrices::Matrix(){}
	Matrix(size_t w, size_t h):shared_matrices::Matrix(w,h){this->clear();}


	/////////////
	// METHODS //
	/////////////

	void mean_row(Matrix& mean) {
		matrix_sumt_float(mean, data, width, height);
		vector_sdiv_float(mean, height, width);
	}

	void centering(Matrix& mu_out, Matrix& out) {
		out = *this;
		mean_row(mu_out);
		for(int i=0; i<height; i++) vector_sub_float(out.get_row(i), mu_out, width);
	}

	void covariance(Matrix& X) {
		if(X.width > 0) {covariance_mt(X); return;}
		this->clear();
		matrix_CpAAt_float(data, X, X.width, X.height);
		(*this) /= X.height;
	}

	void covariance_mt(Matrix& X);


	void eigendecompose(Matrix& U, Matrix& L) {
		DBG_START("eigendecomposition");

		_check_("matrix");

		U = (*this); L.clear();
		matrix_eigSym_float(U, L, width);
		matrix_sortEig_float(U, L, width);

		for(int i=0; i<width; i++) if(isnan(L.data[i])) {
			DBG("ERROR : eigendecompose produced NAN -> recomputing with a little added delta");
			for(int i=0; i<width; i++) data[i*width+i] += 0.00001f;
			eigendecompose(U,L);
			return;
		}

		DBG_END();
	}

	void reconstruct(Matrix& U, Matrix& L) {
		DBG_START("reconstruction");
		this->clear();
		matrix_CpAdAt_float(data, U, L,width,height);
		DBG_END();
	}

	void _check_(const char* what) {
		for(int i=0; i<width; i++) for(int j=0; j<height; j++) {
			if(isnan(data[i*width+j])) FATAL("NAN FAULT : in " << what << " [" << i << "," << j << "]");
		}
	}

	inline float min() { return vector_min_float(data, width*height); }
	inline float max() { return vector_max_float(data, width*height); }

	inline float n2() { return matrix_n2_float(data,width,height); }
	inline float n1() { return vector_n1_float(data,width*height); }
	inline float l2(Matrix& B) { return matrix_l2_float(data,B,height,width); }

	float l2_normalized(Matrix& B) {
		float norm = n2(); if(norm==0) norm=1;
		return l2(B) / norm;
	}

	int count_nonzeros() {
		int n = 0;
		for(int i=0; i<width*height; i++) if(data[i]!=0) n++;
		return n;
	}

	void normalize_l1() { (*this)/=n1(); }


	void AAt(Matrix& A) {	matrix_CpAAt_float(data,A,A.height,A.width); }


	void set_part(int x, int y, Matrix& mat) {
		for(int l=0; l<mat.height; l++)
			memcpy(get_row(l+y)+x, mat.get_row(l), mat.width*sizeof(float));
	}

	void clear_row(int y) {
		memset(get_row(y), 0, sizeof(float)*width);
	}

	void set(int x, int y, float val) {data[x+y*width]=val;}

	void sadd(float s, Matrix& m) {
		for(int x=0; x<width*height; x++) data[x] += s*m[x];
	}

	///////////////
	// OPERATORS //
	///////////////

	void operator=(const int* i) { *((shared_matrices::Matrix*)this) = i; }

	void operator/=(float f) {vector_sdiv_float(data, f, width*height);}
	void operator*=(float f) {vector_smul_float(data, f, width*height);}

	Matrix& operator+=(Matrix & m) { vector_add_float(data, m, width*height); return (*this); }
	Matrix& operator-=(Matrix & m) {
		if(m.height==1) {
			for(int i=0; i<height; i++) vector_sub_float(get_row(i), m, width);
		}
		else vector_sub_float(data, m, width*height);
		return (*this);
	}


	/////////
	// DBG //
	/////////

	void dbg_range() { DBG(min() << " -> " << max()); }
};


#endif /* MATRIX_H_ */
