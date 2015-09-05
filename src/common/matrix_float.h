/*
 * matrix.h
 *
 *  Created on: 15 nov. 2013
 *      Author: jfellus
 */

#ifndef GOSSIP_MATRIX_H_
#define GOSSIP_MATRIX_H_

#include "math.h"

#include "../retin/toolbox/core/SharedMatrix.h"
#include "../retin/toolbox/core/sharedmem/shared_mem.h"
#include "../retin/toolbox/core/sharedmem/matrix_reader.h"
#include "../retin/toolbox/core/sharedmem/matrix_writer.h"
#include "../retin/toolbox/core/string_retin.h"
#include "../retin/toolbox/algebra/matrix_float.h"
#include "../retin/toolbox/algebra/matrix_double.h"
#include "../retin/toolbox/algebra/vector_float.h"
#include "common/utils.h"

using namespace retin;

class MatrixFloat {
public:
	size_t width, height;
	float* data;
	string CREATED_HOW;
public:

	MatrixFloat() {
		data = 0;
		width = height = 0;
		CREATED_HOW = "NOT CREATED !!!";
	}
	MatrixFloat(size_t w, size_t h) {
		data = new float[w*h];
		width = w;
		height = h;
		clear();
		CREATED_HOW = "EMPTY OK";
	}

	MatrixFloat(const MatrixFloat& m) {
	//	DBG("COPY");
		width = m.width;
		height = m.height;
		data = new float[width*height];
		*this = m;
		CREATED_HOW = "COPY";
	}

	virtual ~MatrixFloat() {
	//	DBG("Free : " << CREATED_HOW);
		if(data) delete data;
		data = 0;
	}

	void create_ref(float* data, size_t width, size_t height) {
		this->width = width;
		this->height = height;
		this->data = data;
		CREATED_HOW = "REF";
	}

	void init(size_t w, size_t h) {
		if(data) delete[] data;
		data = new float[w*h];
		width = w;
		height = h;
		CREATED_HOW = "NOT THEN REINITED";
	}


	float* get_row(size_t i) {return data + i*width;}
	operator bool() const { return data!=0; }
	operator float*() const {return data; }
	float& operator[](int i) {return data[i]; }

	/////////////
	// METHODS //
	/////////////

	void clear() {
		if(!data) return;
		memset(data, 0, sizeof(float)*width*height);
	}

	// OK MATLAB
	void identity() {
		clear(); for(int i=0; i<width; i++) data[i*width+i] = 1;
	}

	// OK MATLAB
	void ones() {
		for(int i=0; i<width*height; i++) data[i] = 1;
	}

	// OK MATLAB
	void diag(MatrixFloat& out) {
		for(int i=0; i<MIN(width,height); i++) out.data[i] = data[i*width+i];
	}

	// OK MATLAB
	MatrixFloat diag() {
		if(height==1) {
			MatrixFloat res(width,width); res.clear();
			for(int i=0; i<width; i++) res.data[i*width+i] = data[i];
			return res;
		} else {
			MatrixFloat res(width, 1); diag(res); return res;
		}
	}

	// OK MATLAB
	void mean_row(MatrixFloat& mean) {
		if(!mean.data) { mean.init(width, 1); mean.clear(); }
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				mean[j] += data[i*width + j];
			}
		}
		mean /= height;
	}

	// OK MATLAB
	MatrixFloat mean_row() { MatrixFloat mu; mean_row(mu); return mu; }

	// OK MATLAB
	void substract_meanrow() {
		MatrixFloat mu = mean_row();
		MatrixFloat ones(1,height); ones.ones();
		MatrixFloat z = ones*mu;
		*this -= z;
	}

	// OK MATLAB
	void centering(MatrixFloat& mu_out, MatrixFloat& out) {
		out = *this;
		mean_row(mu_out);
		for(int i=0; i<height; i++) vector_sub_float(out.get_row(i), mu_out, width);
	}

	// OK MATLAB
	void covariance(MatrixFloat& X) {
		if(!data) init(X.width,X.width);
		this->clear();
		matrix_CpAAt_float(data, X, X.width, X.height);
		(*this) /= (X.height);
	}

	// OK MATLAB
	void gram(MatrixFloat& X) {
		if(!data) init(X.height,X.height);
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				data[i*width + j] = 0;
				for(int k=0; k<X.width; k++) {
					data[i*width + j] += X[i*X.width + k]*X[j*X.width + k];
				}
			}
		}
	}

	// OK MATLAB !
	void correlation(MatrixFloat& X) {
		if(!data) init(X.width,X.width);
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				data[i*width + j] = 0;
				for(int k=0; k<X.height; k++) {
					data[i*width + j] += X[k*X.width + i]*X[k*X.width + j];
				}
			}
		}
	}

	// OK MATLAB !
	void QR(MatrixFloat& R);
	void QR(MatrixFloat& Q, MatrixFloat& R);
	void qr(MatrixFloat& Q, MatrixFloat& R);
	void qr(MatrixFloat& R);

	// OK MATLAB !
	void oi(MatrixFloat& U, MatrixFloat& L) {
		if(!U.data) {U.init(width, height);}
		if(!L.data) {L.init(MIN(width,height),1);}
		U.identity();
		MatrixFloat M(U.width, U.height);
		for(int it=0; it<30; it++) {
			M.clear(); M.CpAB(*this, U);
			M.qr(U,L);
		}
	}



	// OK MATLAB = r^(-1/2)
	void m12(MatrixFloat& r) {
		if(!data) init(r.width, r.height);
		for(int i=0; i<r.width*r.height; i++) {
			if(r.data[i]<1e-5) data[i] = 0;
			else data[i] = 1/sqrt(r.data[i]);
		}
	}


	// OK MATLAB !!
	void XUL12(MatrixFloat& X, MatrixFloat& U, MatrixFloat& L) {
		if(!data) init(U.width,X.width);
		for(int i=0; i<MIN(height,X.width); i++) {
			for(int j=0; j<MIN(width,L.width); j++) {
				data[i*width + j] = 0;
				for(int k=0; k<X.height; k++) {
					if(L[j]>1e-5) data[i*width + j] += X[k*X.width + i] * U[k*U.width + j] / sqrt(L[j]);
				}
			}
		}
	}

	void project(MatrixFloat& X, MatrixFloat& U, MatrixFloat& L) {
		if(!data) init(U.width,X.height);
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				data[i*width + j] = 0;
				for(int k=0; k<X.width; k++) data[i*width + j] += X[i*X.width+k]*U[k*U.width +j]*sqrt(L[j]);
			}
		}
	}

	void project(MatrixFloat& X, MatrixFloat& U) {
		if(!data) init(U.width,X.height);
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				data[i*width + j] = 0;
				for(int k=0; k<X.width; k++) data[i*width + j] += X[i*X.width+k]*U[k*U.width +j];
			}
		}
	}

	// OK MATLAB !!
	void unproject(MatrixFloat& X, MatrixFloat& U) {
		if(!data) init(U.height,X.height);
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				data[i*width + j] = 0;
				for(int k=0; k<X.width; k++) data[i*width + j] += X[i*X.width+k]*U[j * U.width + k];
			}
		}
	}

	void normalize_signs() {
		for(int i=0; i<width; i++) {
			if(data[i]>=0) continue;
			for(int j=0; j<height; j++) data[j*width + i] = -data[j*width + i];
		}
	}

	void ULUtQ(MatrixFloat& U, MatrixFloat& L, MatrixFloat& Q) {
		if(!data) {init(Q.width, U.height); clear();}

		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				//data[i*Q.width + j] = 0;
				for(int k=0; k<L.width; k++) {
					for(int l=0; l<Q.height; l++) {
						data[i*width + j] += U[i*U.width + k] * L[k] * U[l*U.width + k] * Q[l*Q.width + j];
					}
				}
			}
		}
	}

	void covariance_mt(MatrixFloat& X);

	void diag_mul_sqrt(MatrixFloat& diag) {
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				data[i*width + j] *= sqrt(diag.data[j]);
			}
		}
	}

	void diag_mulL_sqrt(MatrixFloat& diag) {
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				data[i*width + j] *= sqrt(diag.data[i]);
			}
		}
	}


	void make_h_block(MatrixFloat& A, MatrixFloat& B) {
		if(!data) init(A.width+B.width, A.height);
		for(int i=0; i<height; i++) {
			for(int j=0; j<A.width; j++) data[i*width + j] = A.data[i*A.width + j];
			for(int j=0; j<B.width; j++) data[i*width + A.width + j] = B.data[i*B.width + j];
		}
	}

	void PCA_via_gram(MatrixFloat &U, MatrixFloat& L) {
		MatrixFloat G(height,height);
		MatrixFloat UU(height,height);
		MatrixFloat LL(height, 1);
DDD(		G.gram(*this); )
DDD(		G.eigendecompose(UU, LL); )
		L.clear();
		for(int i=0; i<MIN(L.width,LL.width); i++) L[i] = LL[i];
DDD(		U.XUL12(*this,UU,L); )
	}

	void PCA(MatrixFloat &U, MatrixFloat &L, int q) {
		if(q>width) q = width;
		if(q>height) q = height;
		if(!L.data) {L.init(q,1); L.clear();}
		if(!U.data) {U.init(q,width); U.clear();}

		if(height < width) {
		//	DBG("USE GRAM TRICK");
			PCA_via_gram(U,L);
		}
		else {
			MatrixFloat C(width,width);
			MatrixFloat UU,LL;
	DDD(		C.covariance(*this); )
DDD(			C.eigendecompose(UU,LL); )
			for(int i=0; i<q; i++) {
				L[i] = LL[i];
				for(int j=0; j<width; j++) U[j*U.width+i] = UU[j*UU.width+i];
			}
		}

	}


	// OK MATLAB
	void eigendecompose(MatrixFloat& U, MatrixFloat& L);

	// OK MATLAB
	void transpose(const MatrixFloat& m) {
		if(!data) init(m.height, m.width);
		for(int i=0;i<m.height;i++) {
			for(int j=0; j<m.width; j++) {
				data[j*width+i] = m.data[i*m.width+j];
			}
		}
	}

	// OK MATLAB
	MatrixFloat transpose() const {	MatrixFloat res(height,width); res.transpose(*this); return res;}

	// OK MATLAB
	void reconstruct(MatrixFloat& U, MatrixFloat& L) {
		if(!data) init(U.height, U.height);
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				data[i*width + j] = 0;
				for(int k=0; k<U.width; k++) {
					data[i*width + j] += U[i*U.width + k]*U[j*U.width + k]*L[k];
				}
			}
		}
	}

	void _check_(const char* what) {
		for(int i=0; i<width; i++) for(int j=0; j<height; j++) {
			if(isnan(data[i*width+j])) FATAL("NAN FAULT : in " << what << " [" << i << "," << j << "]");
		}
	}

	inline float min() { return vector_min_float(data, width*height); }
	inline float max() { return vector_max_float(data, width*height); }
	inline float sum() { return vector_sum_float(data, width*height); }

	inline float n2() { return matrix_n2_float(data,width,height); }
	inline float n1() { return vector_n1_float(data,width*height); }
	inline float l2(MatrixFloat& B) { return matrix_l2_float(data,B,height,width); }

	float l2_normalized(MatrixFloat& B) {
		float norm = n2(); if(norm==0) norm=1;
		return l2(B) / norm;
	}

	int count_nonzeros() {
		int n = 0;
		for(int i=0; i<width*height; i++) if(data[i]!=0) n++;
		return n;
	}

	void normalize_l1() { (*this)/=n1(); }


	void AAt(MatrixFloat& A) {	matrix_CpAAt_float(data,A,A.height,A.width); }


	void set_part(int x, int y, MatrixFloat& mat) {
		for(int l=0; l<mat.height; l++)
			memcpy(get_row(l+y)+x, mat.get_row(l), mat.width*sizeof(float));
	}

	void clear_row(int y) {
		memset(get_row(y), 0, sizeof(float)*width);
	}

	void set(int x, int y, float val) {data[x+y*width]=val;}

	void sadd(float s, MatrixFloat& m) {
		for(int x=0; x<width*height; x++) data[x] += s*m[x];
	}

	int nearest_neighbor(float* v) {
		float d,mind=FLT_MAX;
		int nn = -1;
		for(int k=0; k<height; k++) {
			d = vector_l2p2_float(get_row(k),v,width);
			if(d <= mind) {
				mind = d;
				nn = k;
			}
		}
		return nn;
	}

	float nearest_neighbor_distance(float* v) {
		float d,mind=FLT_MAX;
		for(int k=0; k<height; k++) {
			d = vector_l2p2_float(&data[k*width],v,width);
			if(d <= mind) mind = d;
		}
		return mind;
	}

	void row_sdiv(int row, float s) {
		vector_sdiv_float(get_row(row), s, width);
	}

	void row_smul(int row, float s) {
		vector_smul_float(get_row(row), s, width);
	}

	void row_add(int row, float* v) {
		vector_add_float(get_row(row), v, width);
	}


	void CpAB(MatrixFloat& A, MatrixFloat& B) {
		if(!data) {init(B.width, A.height); clear(); }
		for(int i=0; i<height; i++) {
			for(int j=0; j<width; j++) {
				for(int k=0; k<A.width; k++) {
					data[i*width + j] += A.data[i*A.width + k] * B.data[k*B.width + j];
				}
			}
		}
	}


	void CpAtB(MatrixFloat& A, MatrixFloat& B) {
		matrix_CpAtB_float(data, A, B, width, A.height, height);
	}

	void rand(float min = 0, float max = 1);

	void integrate() {
		for(int i=1; i<width*height; i++)
			data[i] = data[i-1] + data[i];
	}


	void randGaussian(int D, int n, int p) {
		MatrixFloat mu(D,1); mu.rand(0,1);
		MatrixFloat sigma(D,D); rand_covariance(sigma,1, p);
		init(D,n);
		for(int i=0; i<n; i++) {
			randvec_gaussian(&data[i*D],mu,sigma);
		}
	}

	///////////////
	// OPERATORS //
	///////////////

	void operator=(const int* i) { *((shared_matrices::Matrix*)this) = i; }
	void operator=(const float* f) { *((shared_matrices::Matrix*)this) = f; }
	float operator=(float x) { *((shared_matrices::Matrix*)this) = x; return x;}
	void operator=(const MatrixFloat& m) {
		if(!data) init(m.width, m.height);
		memcpy(data, m.data, m.width*m.height*sizeof(float));
	}

	void operator/=(float f) {vector_sdiv_float(data, f, width*height);}
	void operator*=(float f) {vector_smul_float(data, f, width*height);}

	MatrixFloat& operator+=(MatrixFloat & m) { vector_add_float(data, m, width*height); return (*this); }
	MatrixFloat& operator-=(MatrixFloat & m) {
		if(m.height==1) {
			for(int i=0; i<height; i++) vector_sub_float(get_row(i), m, width);
		}
		else vector_sub_float(data, m, width*height);
		return (*this);
	}

	MatrixFloat operator+(MatrixFloat &m) {
		MatrixFloat res(width, height); res = *this; res+=m;
		return res;
	}

	MatrixFloat operator-(MatrixFloat &m) {
			MatrixFloat res(width, height); res = *this; res-=m;
			return res;
	}

	MatrixFloat operator*(MatrixFloat &m) {
		MatrixFloat res(m.width, height); res.clear(); res.CpAB(*this, m);
		return res;
	}

	/////////
	// DBG //
	/////////

	void dbg_range() { DBG(min() << " -> " << max()); }

	void dump(size_t nbrows = 0, size_t nbcols = 0) {
		if(nbrows==0) nbrows = height;
		if(nbcols==0) nbcols = width;

		for(size_t i=0; i<nbrows; i++) {
			if(i) putc('\n',stdout);
			for(size_t j=0; j<nbcols; j++) {
				if(j) putc(' ',stdout);
				printf("%f", data[i*width+j]);
			}
			fflush(stdout);
		}
		putc('\n',stdout);
		putc('\n',stdout);
	}

	void dumpdim() {
		std::cout << height << "x" << width << "\n";
	}

	////////
	// IO //
	////////


	inline bool save(const char* file) const {return this->write(file);}
	inline bool save(const string& file) const {return this->write(file.c_str());}
	inline void load(const string& file) {this->read(file.c_str()); }


	bool read(const char* file) {
		matrix_reader<float>* reader = 0;
		std::string file_name = file;
		if (string_has_file_ext(file_name,".fvec")) reader = new fvec_reader(file);
	    else if (string_has_file_ext(file_name,".hvec8")) reader = new hvec8_reader(file);
	    else if (string_has_file_ext(file_name,".csv") || string_has_file_ext(file_name,".m")) reader = new csv_reader<float>(file);
	    else if (string_has_file_ext(file_name,".txt")) reader = new txt_reader<float>(file);
	    else if (string_has_file_ext(file_name,".bin16")) reader = new bin_reader(file);
	    else if (string_has_file_ext(file_name,".peter8")) reader = new peter8_reader(file);
	    else if (string_has_file_ext(file_name,".idx3-ubyte")) reader = new idx3_ubyte_reader(file);

		if(!reader) throw "No reader for this file format";

		height = reader->get_height();
		width = reader->get_width();

		data = new float[width*height];

		fflush(stdout);

	    reader->read(data);

		delete reader;

		return true;
	}

	bool write(const char* file) const {
		matrix_writer<float>* writer = 0;
		std::string file_name = file;
		if (string_has_file_ext(file_name,".fvec")) writer = new fvec_writer(file);
		else if (string_has_file_ext(file_name,".hvec8")) writer = new hvec8_writer(file);
		else if (string_has_file_ext(file_name,".csv") || string_has_file_ext(file_name,".m")) writer = new csv_writer<float>(file);
		else if (string_has_file_ext(file_name,".txt")) writer = new txt_writer<float>(file);
		else if (string_has_file_ext(file_name,".ivecs")) writer = new ivecs_writer(file);
		else if (string_has_file_ext(file_name,".peter8")) writer = new peter8_writer(file);

		if(!writer) throw "No writer for this file format";

		writer->write(data, width, height);
		delete writer;
		return true;
	}

};




#endif /* MATRIX_H_ */
