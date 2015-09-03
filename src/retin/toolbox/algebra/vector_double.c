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
 * \file vector.c
 * \author Philippe H. Gosselin
 * \version 1.0
 */

#include "vector_double.h"

double vector_n1_double(const double* v, size_t n) {
	size_t i;
	double s, x;
	s = 0;
	for (i = 0; i < n; i++) {
		x = v[i];
		if (x > 0)
			s += v[i];
		else
			s -= v[i];
	}
	return s;
}

double vector_n2_double(const double* v, size_t n) {
	return sqrt(vector_n2p2_double(v, n));
}

double vector_n2p2_double(const double* v, size_t n) {
	return vector_ps_double(v, v, n);
}

double vector_l2_double(const double* v1, const double* v2, size_t n) {
	return sqrt(vector_l2p2_double(v1, v2, n));
}

double vector_l2p2_double(const double* v1, const double* v2, size_t n) {
	size_t i;
	double s, x;
	s = 0;
	for (i = 0; i < n; i++) {
		x = v1[i] - v2[i];
		s += x * x;
	}
	return s;
}

#if defined(ALGEBRA_AVX)
ALGEBRA_INLINE double vector_ps_double (const double* pa,const double* pb,size_t n) {
	if(ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
		size_t q = n/4;
		size_t r = n%4;
		double w = 0;

		if(q>0) {
			__m256d acc = _mm256_setzero_pd();
			__m256d i1 = _mm256_load_pd(pa);
			__m256d j1 = _mm256_load_pd(pb);
			pa += 4;
			pb += 4;
			__m256d s = _mm256_mul_pd(i1, j1);
			acc = _mm256_add_pd(acc, s);

			while(--q != 0) {
				// load
				i1 = _mm256_load_pd(pa);
				j1 = _mm256_load_pd(pb);
				pa += 4;
				pb += 4;
				// multiplie
				s = _mm256_mul_pd(i1, j1);
				// accumule
				acc = _mm256_add_pd(acc, s);
			}
			// sum finale
			// add horizontal
			acc = _mm256_hadd_pd(acc, acc);
			// échange 128bits haut et bas
			__m256d accp = _mm256_permute2f128_pd(acc, acc, 1);
			// add vertical
			acc = _mm256_add_pd(acc, accp);
			// extract
			_mm_store_sd(&w, _mm256_extractf128_pd(acc,0));
		}
		return w + vector_ps_double_basic(pa, pb, r);
	}
	return vector_ps_double_basic(pa, pb, n);
}

#elif defined(ALGEBRA_SSSE3)

double vector_ps_double (const double* pa,const double* pb,size_t n)
{
	size_t k;
	/* multiplication 4 par 4 */
	size_t q = n / 4;
	size_t r = n % 4;
	double w;
	_mm_prefetch (pa,_MM_HINT_NTA);
	_mm_prefetch (pb,_MM_HINT_NTA);
	if (q > 0) {
		__m128d acc1 = _mm_setzero_pd();
		__m128d acc2 = _mm_setzero_pd();
		if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
			for (k=0;k<q;k++) {
				/* Charge 2 doubles dans chaque tableau */
				__m128d i1 = _mm_load_pd(pa);
				__m128d j1 = _mm_load_pd(pb);
				__m128d i2 = _mm_load_pd(pa+2);
				__m128d j2 = _mm_load_pd(pb+2);
				/* incrément de 4 doubles en tout (2 pour i et 2 pour j) */
				/* Multiplie */
				__m128d s1 = _mm_mul_pd(i1,j1);
				__m128d s2 = _mm_mul_pd(i2,j2);
				pa += 4;
				pb += 4;
				/* Accumule */
				acc1 = _mm_add_pd(acc1,s1);
				acc2 = _mm_add_pd(acc2,s2);
			}
		}
		else {
			for (k=0;k<q;k++) {
				/* Charge 2 doubles dans chaque tableau */
				__m128d i1 = _mm_loadu_pd(pa);
				__m128d j1 = _mm_loadu_pd(pb);
				__m128d i2 = _mm_loadu_pd(pa+2);
				__m128d j2 = _mm_loadu_pd(pb+2);
				/* Multiplie */
				__m128d s1 = _mm_mul_pd(i1,j1);
				__m128d s2 = _mm_mul_pd(i2,j2);
				pa += 4;
				pb += 4;
				/* Accumule */
				acc1 = _mm_add_pd(acc1,s1);
				acc2 = _mm_add_pd(acc2,s2);
			}
		}
		/* Somme finale */
		acc1 = _mm_add_pd(acc1,acc2);
		acc1 = _mm_hadd_pd(acc1,acc1);
		_mm_store_sd(&w,acc1);
	}
	else {
		w = 0;
	}
	for (k=0;k<r;k++)
	w += (*pa++) * (*pb++);
	return w;
}

#else

double vector_ps_double(const double* v1, const double* v2, size_t n) {
	return vector_ps_double_basic(v1, v2, n);
}

#endif

double vector_ps_double_basic(const double* v1, const double* v2, size_t n) {
	size_t i;
	double s = 0;
	for (i = 0; i < n; i++)
		s += v1[i] * v2[i];
	return s;
}

double vector_max_double(const double* v, size_t n) {
	size_t i;
	double x = v[0];
	for (i = 1; i < n; i++)
		if (x < v[i])
			x = v[i];
	return x;
}

void vector_linear_double(double* v, const double* v1, const double alpha,
		const double* v2, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v[i] = v1[i] + alpha * v2[i];
}

void vector_add_stats_double(double* mean, double* var, double lambda,
		const float* v, size_t n) {
	size_t i;
	for (i = 0; i < n; i++) {
		double x = v[i];
		mean[i] += lambda * x;
		var[i] += lambda * x * x;
	}
}

void vector_scpy_double(double* v, double r, size_t n) {
	if (r == 0) {
		memset(v, 0, n * sizeof(double));
	} else {
		size_t i;
		for (i = 0; i < n; i++)
			v[i] = r;
	}
}

void vector_sadd_double(double* v, double r, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v[i] += r;
}

void vector_ssub_double(double* v, double r, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v[i] -= r;
}

void vector_smul_double(double* v, double r, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v[i] *= r;
}

void vector_sdiv_double(double* v, double r, size_t n) {
	size_t i;
	if (r == 0)
		return;
	r = 1.0 / r;
	for (i = 0; i < n; i++)
		v[i] *= r;
}

void vector_add_double(double* v1, const double* v2, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v1[i] += v2[i];
}

#if defined(ALGEBRA_AVX)

ALGEBRA_INLINE void vector_addm_double_aligned_32 (double* v1,double lambda,const double* v2,size_t n)
{
	size_t k;

	__m256d l1 = _mm256_broadcast_sd(&lambda);
	__m256d l2 = _mm256_broadcast_sd(&lambda);
	__m256d l3 = _mm256_broadcast_sd(&lambda);
	__m256d l4 = _mm256_broadcast_sd(&lambda);

	size_t q = n / 16;
	size_t r = n % 16;
	if(q > 0) {
		if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
			for (k=0;k<q;k++) {
				/* Charge 4 valeurs de chaque tableau */
				__m256d i1 = _mm256_load_pd(v1);
				__m256d j1 = _mm256_load_pd(v2);
				__m256d i2 = _mm256_load_pd(v1+4);
				__m256d j2 = _mm256_load_pd(v2+4);
				__m256d i3 = _mm256_load_pd(v1+8);
				__m256d j3 = _mm256_load_pd(v2+8);
				__m256d i4 = _mm256_load_pd(v1+12);
				__m256d j4 = _mm256_load_pd(v2+12);
				/* multiplie */
				j1 = _mm256_mul_pd(j1, l1);
				j2 = _mm256_mul_pd(j2, l2);
				j3 = _mm256_mul_pd(j3, l3);
				j4 = _mm256_mul_pd(j4, l4);
				/* Additionne */
				i1 = _mm256_add_pd(i1,j1);
				i2 = _mm256_add_pd(i2,j2);
				i3 = _mm256_add_pd(i3,j3);
				i4 = _mm256_add_pd(i4,j4);
				/* Sauvegarde */
				_mm256_store_pd(v1, i1);
				_mm256_store_pd(v1+4, i2);
				_mm256_store_pd(v1+8, i3);
				_mm256_store_pd(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
		else {
			for (k=0;k<q;k++) {
				/* Charge 4 valeurs de chaque tableau */
				__m256d i1 = _mm256_loadu_pd(v1);
				__m256d j1 = _mm256_loadu_pd(v2);
				__m256d i2 = _mm256_loadu_pd(v1+4);
				__m256d j2 = _mm256_loadu_pd(v2+4);
				__m256d i3 = _mm256_loadu_pd(v1+8);
				__m256d j3 = _mm256_loadu_pd(v2+8);
				__m256d i4 = _mm256_loadu_pd(v1+12);
				__m256d j4 = _mm256_loadu_pd(v2+12);
				/* multiplie */
				j1 = _mm256_mul_pd(j1, l1);
				j2 = _mm256_mul_pd(j2, l2);
				j3 = _mm256_mul_pd(j3, l3);
				j4 = _mm256_mul_pd(j4, l4);
				/* Additionne */
				i1 = _mm256_add_pd(i1,j1);
				i2 = _mm256_add_pd(i2,j2);
				i3 = _mm256_add_pd(i3,j3);
				i4 = _mm256_add_pd(i4,j4);
				/* Sauvegarde */
				_mm256_storeu_pd(v1, i1);
				_mm256_storeu_pd(v1+4, i2);
				_mm256_storeu_pd(v1+8, i3);
				_mm256_storeu_pd(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
	}

	for(k = 0; k<r; k++)
	v1[k] += lambda*v2[k];
}

#elif defined(ALGEBRA_SSSE3)

ALGEBRA_INLINE void vector_addm_double_aligned_32 (double* v1,double lambda,const double* v2,size_t n)
{
	size_t k;

	__m128d l1 = _mm_load1_pd(&lambda);

	size_t q = n / 2;
	size_t r = n % 2;
	if(q > 0) {
		if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
			for (k=0;k<q;k++) {
				/* Charge 2 valeurs de chaque tableau */
				__m128d i1 = _mm_load_pd(v1);
				__m128d j1 = _mm_load_pd(v2);
				/* multiplie */
				j1 = _mm_mul_pd(j1, l1);
				/* additionne */
				i1 = _mm_add_pd(i1,j1);
				/* Sauvegarde */
				_mm_store_pd(v1, i1);
				v1 += 2;
				v2 += 2;
			}
		}
		else {
			for (k=0;k<q;k++) {
				/* Charge 8 valeurs de chaque tableau */
				__m128d i1 = _mm_loadu_pd(v1);
				__m128d j1 = _mm_loadu_pd(v2);
				j1 = _mm_mul_pd(j1, l1);
				/* Soustrait */
				i1 = _mm_add_pd(i1,j1);
				/* Sauvegarde */
				_mm_storeu_pd(v1, i1);
				v1 += 2;
				v2 += 2;
			}
		}
	}

	for(k = 0; k<r; k++)
	v1[k] += lambda*v2[k];

}

#else

ALGEBRA_INLINE void vector_addm_double_aligned_32(double* v1, double lambda,
		const double* v2, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v1[i] += lambda * v2[i];
}

#endif

void vector_addm_double(double* v1, double lambda, const double* v2, size_t n) {
	vector_addm_double_aligned_32(v1, lambda, v2, n);
}

void vector_addm_double_basic(double* v1, double lambda, const double* v2,
		size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v1[i] += lambda * v2[i];
}

void vector_sub_double(double* v1, const double* v2, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v1[i] -= v2[i];
}

void vector_mul_double(double* v1, double* v2, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v1[i] *= v2[i];
}

void vector_div_double(double* v1, double* v2, size_t n) {
	size_t i;
	for (i = 0; i < n; i++)
		v1[i] /= v2[i];
}

void vector_exp_proba_double(double* p, size_t n) {
	size_t i;
	double x, safe_sum = 0;
	double p_max = vector_max_double(p, n);
	// log( sum_j exp(gamma[j]) ) 
	// = log( sum_j exp(gamma[j] - gamma_max + gamma_max) )
	// = log( sum_j exp(gamma[j] - gamma_max ) * exp(gamma_max) )
	// = log( sum_j exp(gamma[j] - gamma_max ) + log( exp(gamma_max) )
	// = log( sum_j exp(gamma[j] - gamma_max ) + gamma_max
	for (i = 0; i < n; i++) {
		x = exp(p[i] - p_max); // on minimize la valeur à passer dans l'exponentielle
		safe_sum += x;
		p[i] = x;
	}
	vector_sdiv_double(p, safe_sum, n);
}

double vector_min_double(const double* v, size_t n) {
	size_t i;
	double x = v[0];
	for (i = 1; i < n; i++)
		if (x > v[i])
			x = v[i];
	return x;
}

double vector_sum_double(const double* v, size_t n) {
	double w = 0;
	size_t i;
	for (i = 0; i < n; i++)
		w += v[i];
	return w;
}

