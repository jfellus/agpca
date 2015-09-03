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
 * \file vector.h
 * \author Philippe H. Gosselin
 * \version 1.0
 */

#ifndef __algebra_vector_hpp__
#define __algebra_vector_hpp__

#include "vector_uchar.h"
#include "vector_short.h"
#include "vector_float.h"
#include "vector_double.h"

#include "matrix_float.h"

template<typename T> void vector_zero(T* v, size_t n) {
	memset(v, 0, n * sizeof(T));
}
template<typename T> void vector_cpy(T* v1, const T* v2, size_t n) {
	memcpy(v1, v2, n * sizeof(T));
}
template<typename U, typename T> void vector_convert(U* v1, const T* v2,
		size_t n);

template<typename U, typename T> U vector_ps(const T* v1, const T* v2,
		size_t n) {
	throw "Not implemented";
}
template<> inline unsigned int vector_ps(const unsigned char* v1,
		const unsigned char* v2, size_t n) {
	return vector_ps_uchar(v1, v2, n);
}
template<> inline float vector_ps(const unsigned char* v1,
		const unsigned char* v2, size_t n) {
	return vector_ps_uchar(v1, v2, n);
}
template<> inline double vector_ps(const unsigned char* v1,
		const unsigned char* v2, size_t n) {
	return vector_ps_uchar(v1, v2, n);
}
template<> inline float vector_ps(const float* v1, const float* v2, size_t n) {
	return vector_ps_float(v1, v2, n);
}
template<> inline double vector_ps(const float* v1, const float* v2, size_t n) {
	return matrix_ps_float(v1, v2, n, 1);
}

template<typename T> void vector_pow(T* v, float p, size_t n) {
	throw "Not implemented";
}
template<> inline void vector_pow(float* v, float p, size_t n) {
	vector_pow_float(v, p, n);
}

template<typename T> void vector_sdiv(T* v, double r, size_t n) {
	throw "Not implemented";
}
template<> inline void vector_sdiv(float* v, double r, size_t n) {
	vector_sdiv_float(v, r, n);
}
template<> inline void vector_sdiv(double* v, double r, size_t n) {
	vector_sdiv_double(v, r, n);
}

//////////////////////////////////////////////////////////////////////////
template<typename U, typename T>
void vector_convert(U* v1, const T* v2, size_t n) {
	for (size_t i = 0; i < n; i++) {
		v1[i] = (U) v2[i];
	}
}

#endif
