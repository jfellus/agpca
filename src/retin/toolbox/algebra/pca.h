/*
 Copyright © CNRS 2014. 
 Authors: David Picard, Philippe-Henri Gosselin, Romain Negrel, Hedi 
 Tabia, Jerome Fellus
 Contact: picard@ensea.fr

 This software is governed by the CeCILL license under French law and
 abiding by the rules of distribution of free software.  You can  use, 
 modify and/ or redistribute the software under the terms of the CeCILL
 license as circulated by CEA, CNRS and INRIA at the following URL
 "http://www.cecill.info". 

 As a counterpart to the access to the source code and rights to copy,
 modify and redistribute granted by the license, users are provided only
 with a limited warranty  and the software's author,  the holder of the
 economic rights,  and the successive licensors  have only  limited
 liability. 

 In this respect, the user's attention is drawn to the risks associated
 with loading,  using,  modifying and/or developing or reproducing the
 software by the user in light of its specific status of free software,
 that may mean  that it is complicated to manipulate,  and  that  also
 therefore means  that it is reserved for developers  and  experienced
 professionals having in-depth computer knowledge. Users are therefore
 encouraged to load and test the software's suitability as regards their
 requirements in conditions enabling the security of their systems and/or 
 data to be ensured and,  more generally, to use and operate it in the 
 same conditions as regards security. 

 The fact that you are presently reading this means that you have had
 knowledge of the CeCILL license and that you accept its terms.

 */
/**
 * \file PCA.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __algebra_PCA_h__
#define __algebra_PCA_h__

#include "vector_double.h"
#include "matrix_double.h"

#include <iostream>
#include <stdexcept>
#include <vector>

#ifdef RETIN_ENABLE_BOOST
#include <boost/thread.hpp>
#endif

namespace algebra {

class PCA {

public:
	class Eigen {
		double l;
		std::vector<double> v;
	public:
		Eigen(double l = 0, const double* pv = NULL, size_t n = 0) :
				l(l), v(n) {
			for (size_t i = 0; i < n; i++) {
				v[i] = pv[i];
			}
		}
		bool operator<(const Eigen& e) const {
			return fabs(l) > fabs(e.l);
		}
		double getValue() const {
			return l;
		}
		const double* getVector() const {
			return &v[0];
		}
	};

	PCA(bool verbose = false) :
			verbose(verbose), dim(0), count(0) {
	}

	void setVerbose(bool b) {
		verbose = b;
	}

	void init(size_t d) {
		count = 0;
		dim = d;
		temp.clear();
		temp.resize(dim);
		mean.clear();
		mean.resize(dim);
		covs.clear();
		covs.resize(dim * dim);
	}

	void add(const double* x, size_t xdim) {
		if (dim == 0)
			init(xdim);
		if (xdim != dim)
			throw std::runtime_error("PCA::add() Invalid dim");
		if (verbose && (count % 1000) == 0)
			std::cout << "..." << count << "..." << std::endl;
		vector_add_double(&mean[0], x, dim);
		matrix_Cpaat_double(&covs[0], x, dim);
		count++;
	}

	void add(const float* x, size_t xdim) {
		if (dim == 0)
			init(xdim);
		if (xdim != dim)
			throw std::runtime_error("PCA::add() Invalid dim");
		for (size_t i = 0; i < xdim; i++)
			temp[i] = x[i];
		add(&temp[0], xdim);
	}

	void center();
	void run(size_t maxDim = 0, double maxEnergy = 100);

	size_t getDim() const {
		return dim;
	}
	size_t getCount() const {
		return count;
	}
	double getEnergy() const {
		return totEnergy;
	}
	const double* getMean() const {
		return &mean[0];
	}
	const double* getCovariances() const {
		return &covs[0];
	}
	const std::vector<Eigen>& getEigens() const {
		return eigens;
	}

	//! Calcule toutes les valeurs propres. A utiliser à la place de "run".
	void computeEigens();

private:
	bool verbose;
	size_t dim, count;
	double totEnergy;
	std::vector<double> temp, mean, covs;
	std::vector<Eigen> eigens;

	void sort(size_t maxDim, double maxEnergy = 100);

public:
	void project(float* output, const float* input, const float* mean,
			size_t dim, const float* eigenvectors, const float* eigenvalues,
			size_t eigencount) {
		/*int dictSize = dict->size();
		 int featureDim = dict->dim();
		 vector<float> temp()
		 vector_sub_float(input->data(), mean->data(), featureDim);
		 for (int k = 0; k < dictSize; k++) {
		 float w = 1 / sqrt((*vars)[k]);
		 (*output)(k) = w * vector_ps_float(input->data(), dict->data(k), featureDim);
		 }*/
	}
	static void matrixPca(const double *cov, double *vectors, double *values,
			size_t dim) {
		PCA pca;
		pca.init(dim);
		memcpy(&(pca.covs[0]), cov, dim * dim * sizeof(double));
		pca.count = 1;
		pca.run();
		memset(vectors, 0, dim * dim * sizeof(double));
		memset(values, 0, dim * sizeof(double));
		for (size_t i = 0; i < pca.eigens.size(); i++) {
			values[i] = pca.eigens[i].getValue();
			memcpy(vectors + dim * i, pca.eigens[i].getVector(),
					dim * sizeof(double));
		}
	}
};

#ifdef RETIN_ENABLE_BOOST
class SynchronizedPCA {
	PCA pca;
	mutable boost::mutex the_mutex;
public:
	void setVerbose(bool b) {
		boost::mutex::scoped_lock lock(the_mutex);
		pca.setVerbose(b);
	}
	void add (const double* x,size_t xdim) {
		boost::mutex::scoped_lock lock(the_mutex);
		pca.add(x,xdim);
	}
	void center() {
		boost::mutex::scoped_lock lock(the_mutex);
		pca.center();
	}
	void run(int maxDim,double maxEnergy) {
		boost::mutex::scoped_lock lock(the_mutex);
		pca.run(maxDim,maxEnergy);
	}
	size_t getDim() const {return pca.getDim();}
	size_t getCount() const {return pca.getCount();}
	const double* getMean() const {return pca.getMean();}
	const double* getCovariances() const {return pca.getCovariances();}
	const std::vector<PCA::Eigen>& getEigens() const {return pca.getEigens();}
};
#endif

}

#endif
