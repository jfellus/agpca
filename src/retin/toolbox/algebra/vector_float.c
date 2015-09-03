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

#include "vector_float.h"

void    vector_rand_float (float* v,float fMin,float fMax,size_t n)
{
    size_t i;
    double r;
    for (i=0;i<n;i++) {
        r = rand()%1000000000;
        r /= 1000000000;
        v[i] = fMin + r * (fMax - fMin);
    }
}


void	vector_linear_float (float* v,const float* v1,const float alpha,const float* v2,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
	    v[i] = v1[i]+alpha*v2[i];
}

void	vector_scpy_float (float* v,float r,size_t n)
{
	if (r == 0)
	{
		memset(v,0,n*sizeof(float));
	}
	else
	{
		size_t i;
		for (i=0;i<n;i++)
			v[i] = r;
	}
}

void	vector_sadd_float (float* v,float r,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v[i] += r;
}

void	vector_ssub_float (float* v,float r,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v[i] -= r;
}

void	vector_smul_float (float* v,float r,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v[i] *= r;
}

void	vector_sdiv_float (float* v,float r,size_t n)
{
	size_t i;
    if (r == 0)
        return;
    r = 1.0f / r;
	for (i=0;i<n;i++)
		v[i] *= r;
}

float	vector_max_float (float* v,size_t n)
{
	size_t i;
	float x = v[0];
	for (i=1;i<n;i++)
		if (x < v[i]) x = v[i];
	return x;
}

float	vector_min_float (float* v,size_t n)
{
	size_t i;
	float x = v[0];
	for (i=1;i<n;i++)
		if (x > v[i]) x = v[i];
	return x;
}

size_t	vector_argmax_float (float* v,size_t n)
{
	size_t i,j = 0;
	float x = v[0];
	for (i=1;i<n;i++)
		if (x < v[i]) { x = v[i]; j = i; }
	return j;
}

size_t	vector_argmin_float (float* v,size_t n)
{
	size_t i,j = 0;
	float x = v[0];
	for (i=1;i<n;i++)
		if (x > v[i]) { x = v[i]; j = i; }
	return j;
}
 
void	vector_rescale_float (float* v,float rMin,float rMax,size_t n)
{
	size_t i;
	float vMin = v[0];
	float vMax = v[0];
	for (i=1;i<n;i++)
	{
		if (vMin > v[i]) vMin = v[i];
		if (vMax < v[i]) vMax = v[i];
	}
	for (i=0;i<n;i++)
		v[i] = (rMax - rMin)*(v[i] - vMin)/(vMax - vMin) + rMin;
}

void	vector_abs_float (float* v,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		if (v[i] < 0) v[i] = -v[i];
}

#if defined(ALGEBRA_SSE2)

void	vector_sqrt_float (float* pv,size_t n)
{
    size_t k;
    size_t q = n / 4;
    size_t r = n % 4;
    static const unsigned int mask_abs = 0x7FFFFFFF;
    const __m128 _mm_mask_abs_ps = _mm_load1_ps((const float*)&mask_abs);
    static const unsigned int mask_neg = 0x80000000;
    const __m128 _mm_mask_neg_ps = _mm_load1_ps((const float*)&mask_neg);
    const __m128 zero = _mm_setzero_ps();
    if (q > 0) {
	if (ALGEBRA_IS_ALIGNED(pv)) {
	    for (k=0;k<q;k++) {
		/* Charge valeurs */
		__m128 i1 = _mm_load_ps(pv);
                /* Valeur absolue  */
		__m128 a1 = _mm_and_ps(i1,_mm_mask_abs_ps);
		/* Racine carrée */
		__m128 r1 = _mm_sqrt_ps(a1);
                /* Entrée négative ? */
                __m128 s1 = _mm_cmplt_ps(i1,zero); /* s1 = 0xFFFFF si négatif */
                s1 = _mm_and_ps(s1,_mm_mask_neg_ps); /* s1 = 0x8000000 si négatif */
                r1 = _mm_or_ps(r1,s1); /* met le signe négatif si s1 négatif */
                /* Sauvegarde */
                _mm_store_ps(pv,r1);
		pv += 4;
	    }
	}
	else {
	    for (k=0;k<q;k++) {
		/* Charge valeurs */
		__m128 i1 = _mm_loadu_ps(pv);
                /* Valeur absolue */
		__m128 a1 = _mm_and_ps(i1,_mm_mask_abs_ps);
		/* Racine carrée */
		__m128 r1 = _mm_sqrt_ps(a1);
                /* Entrée négative ? */
                __m128 s1 = _mm_cmplt_ps(i1,zero); /* s1 = 0xFFFFF si négatif */
                s1 = _mm_and_ps(s1,_mm_mask_neg_ps); /* s1 = 0x8000000 si négatif */
                r1 = _mm_or_ps(r1,s1); /* met le signe négatif si s1 négatif */
                /* Sauvegarde */
                _mm_storeu_ps(pv,r1);
		pv += 4;
	    }
	}
    }
    for (k=0;k<r;k++) {
        float x = *pv;
        if (x >= 0)
            *pv++ = sqrt(x);
        else
            *pv++ = -sqrt(-x);
    }
}

#else

void	vector_sqrt_float (float* v,size_t n)
{
    size_t i;
    for (i=0;i<n;i++) {
        float x = v[i];
        if (x >= 0)
            v[i] = sqrt(x);
        else
            v[i] = -sqrt(-x);
    }
}

#endif

void	vector_pow_float (float* v,float p,size_t n)
{
	size_t i;
        if (p == 1) {
            return;
        }
        if (p == 0.5) {
            vector_sqrt_float(v,n);
            return;
        }
        for (i=0;i<n;i++) {
                float x = v[i];
                if (x >= 0)
                        v[i] = powf(x, p);
                else
                        v[i] = -powf(-x, p);
        }
}


#if defined(ALGEBRA_AVX)

ALGEBRA_INLINE float vector_n1_float_aligned_32(const float * v, size_t q) {
    float n;
    __m256 acc = _mm256_setzero_ps();

    __m256 xplus = _mm256_load_ps(v);
    __m256 xmins = _mm256_sub_ps(_mm256_setzero_ps(), xplus);
    __m256 xabs  = _mm256_max_ps(xplus, xmins);
    acc = _mm256_add_ps(acc, xabs);
    v += 8;
    while(--q!=0) {
        xplus = _mm256_load_ps(v);
        xmins = _mm256_sub_ps(_mm256_setzero_ps(), xplus);
        xabs  = _mm256_max_ps(xplus, xmins);
        acc = _mm256_add_ps(acc, xabs);
        v += 8;
    }
    
    __m256 accp = _mm256_permute2f128_ps(acc, acc, 1);
    acc = _mm256_add_ps(acc, accp);
    acc = _mm256_hadd_ps(acc, acc);
    acc = _mm256_hadd_ps(acc, acc);
    _mm_store_ss(&n,  _mm256_extractf128_ps(acc,0));
    return n;
}

#else
float vector_n1_float_aligned_32(const float * v, size_t q) {
    return vector_n1_float_basic(v, q);
}

#endif

float	vector_n1_float_basic (const float* v,size_t n)
{
	size_t i;
	float s,x;
	s = 0;
	for (i=0;i<n;i++)
	{
		x = v[i];
		if (x > 0) s += v[i];
		else	   s -= v[i];
	}
	return s;
}


float	vector_n1_float (const float * v, size_t n) {
    if (ALGEBRA_IS_ALIGNED(v)) {
        size_t q = n / 8;
        size_t r = n % 8;
        float w = 0;
        if (q > 0) w = vector_n1_float_aligned_32 (v,q);
        return w + vector_n1_float_basic (v+q*8,r);
    }
    return vector_n1_float_basic(v,n);
}

float	vector_n2p2_float_basic (const float* v,size_t n)
{
	size_t i;
	float x;
        double s = 0;
	for (i=0;i<n;i++) {
            x = v[i];
            s += x*x;
        }
	return s;
}

float vector_n2p2_float(const float * v, size_t n) {
    return vector_ps_float(v, v, n);
}


float	vector_n2_float (const float* v,size_t n)
{
    return sqrt(vector_n2p2_float(v,n));
}
   
ALGEBRA_INLINE float	vector_ps_float_basic (const float* v1,const float* v2,size_t n)
{
	size_t i;
	float s = 0;
	for (i=0;i<n;i++)
		s += v1[i]*v2[i];
	return s;
}

#if defined(ALGEBRA_AVX)
 
ALGEBRA_INLINE float	vector_ps_float_aligned_32 (const float* pa,const float* pb,size_t q)
{
    float w;
    __m128 lo,hi;
    __m256 acc1 = _mm256_mul_ps(_mm256_load_ps(pa),_mm256_load_ps(pb));
    __m256 acc2 = _mm256_mul_ps(_mm256_load_ps(pa+8),_mm256_load_ps(pb+8));
    __m256 acc3 = _mm256_mul_ps(_mm256_load_ps(pa+16),_mm256_load_ps(pb+16));
    __m256 acc4 = _mm256_mul_ps(_mm256_load_ps(pa+24),_mm256_load_ps(pb+24));
    pa += 32;
    pb += 32;
    while(--q!=0) {
        __m256 i1 = _mm256_mul_ps(_mm256_load_ps(pa),_mm256_load_ps(pb));
        __m256 i2 = _mm256_mul_ps(_mm256_load_ps(pa+8),_mm256_load_ps(pb+8));
        __m256 i3 = _mm256_mul_ps(_mm256_load_ps(pa+16),_mm256_load_ps(pb+16));
        __m256 i4 = _mm256_mul_ps(_mm256_load_ps(pa+24),_mm256_load_ps(pb+24));
        pa += 32;
        pb += 32;
        acc1 = _mm256_add_ps(acc1,i1);
        acc2 = _mm256_add_ps(acc2,i2);
        acc3 = _mm256_add_ps(acc3,i3);
        acc4 = _mm256_add_ps(acc4,i4);
    }
    acc1 = _mm256_add_ps(acc1,acc2);    
    acc3 = _mm256_add_ps(acc3,acc4);    
    acc1 = _mm256_hadd_ps(acc1,acc3);
    acc1 = _mm256_hadd_ps(acc1,acc1);    
    acc1 = _mm256_hadd_ps(acc1,acc1);    
    lo = _mm256_extractf128_ps (acc1,0);
    hi = _mm256_extractf128_ps (acc1,1);
    _mm_store_ss(&w, _mm_add_ss(lo,hi) );
    return w;
}

ALGEBRA_INLINE float	vector_ps_float_aligned_64 (const float* pa,const float* pb,size_t q)
{
    float w;
    __m128 lo,hi;
    __m256 acc1 = _mm256_mul_ps(_mm256_load_ps(pa),_mm256_load_ps(pb));
    __m256 acc2 = _mm256_mul_ps(_mm256_load_ps(pa+8),_mm256_load_ps(pb+8));
    __m256 acc3 = _mm256_mul_ps(_mm256_load_ps(pa+16),_mm256_load_ps(pb+16));
    __m256 acc4 = _mm256_mul_ps(_mm256_load_ps(pa+24),_mm256_load_ps(pb+24));
    __m256 acc5 = _mm256_mul_ps(_mm256_load_ps(pa+32),_mm256_load_ps(pb+32));
    __m256 acc6 = _mm256_mul_ps(_mm256_load_ps(pa+40),_mm256_load_ps(pb+40));
    __m256 acc7 = _mm256_mul_ps(_mm256_load_ps(pa+48),_mm256_load_ps(pb+48));
    __m256 acc8 = _mm256_mul_ps(_mm256_load_ps(pa+56),_mm256_load_ps(pb+56));
    pa += 64;
    pb += 64;
    while(--q!=0) {
        __m256 i1 = _mm256_mul_ps(_mm256_load_ps(pa),_mm256_load_ps(pb));
        __m256 i2 = _mm256_mul_ps(_mm256_load_ps(pa+8),_mm256_load_ps(pb+8));
        __m256 i3 = _mm256_mul_ps(_mm256_load_ps(pa+16),_mm256_load_ps(pb+16));
        __m256 i4 = _mm256_mul_ps(_mm256_load_ps(pa+24),_mm256_load_ps(pb+24));
        __m256 i5 = _mm256_mul_ps(_mm256_load_ps(pa+32),_mm256_load_ps(pb+32));
        __m256 i6 = _mm256_mul_ps(_mm256_load_ps(pa+40),_mm256_load_ps(pb+40));
        __m256 i7 = _mm256_mul_ps(_mm256_load_ps(pa+48),_mm256_load_ps(pb+48));
        __m256 i8 = _mm256_mul_ps(_mm256_load_ps(pa+56),_mm256_load_ps(pb+56));
        pa += 64;
        pb += 64;
        acc1 = _mm256_add_ps(acc1,i1);
        acc2 = _mm256_add_ps(acc2,i2);
        acc3 = _mm256_add_ps(acc3,i3);
        acc4 = _mm256_add_ps(acc4,i4);
        acc5 = _mm256_add_ps(acc5,i5);
        acc6 = _mm256_add_ps(acc6,i6);
        acc7 = _mm256_add_ps(acc7,i7);
        acc8 = _mm256_add_ps(acc8,i8);
    }
    acc1 = _mm256_add_ps(acc1,acc2);    
    acc3 = _mm256_add_ps(acc3,acc4);    
    acc5 = _mm256_add_ps(acc5,acc6);    
    acc7 = _mm256_add_ps(acc7,acc8);    
    acc1 = _mm256_add_ps(acc1,acc3);    
    acc5 = _mm256_add_ps(acc5,acc7);    
    acc1 = _mm256_add_ps(acc1,acc5);  
    acc1 = _mm256_hadd_ps(acc1,acc1);    
    acc1 = _mm256_hadd_ps(acc1,acc1);    
    lo = _mm256_extractf128_ps (acc1,0);
    hi = _mm256_extractf128_ps (acc1,1);
    _mm_store_ss(&w, _mm_add_ss(lo,hi) );
    return w;
}

#elif defined(ALGEBRA_SSE2)
 
ALGEBRA_INLINE float	vector_ps_float_aligned_32 (const float* pa,const float* pb,size_t q)
{
    float w;
    __m128 acc1 = _mm_mul_ps(_mm_load_ps(pa),_mm_load_ps(pb));
    __m128 acc2 = _mm_mul_ps(_mm_load_ps(pa+4),_mm_load_ps(pb+4));
    __m128 acc3 = _mm_mul_ps(_mm_load_ps(pa+8),_mm_load_ps(pb+8));
    __m128 acc4 = _mm_mul_ps(_mm_load_ps(pa+12),_mm_load_ps(pb+12));
    __m128 acc5 = _mm_mul_ps(_mm_load_ps(pa+16),_mm_load_ps(pb+16));
    __m128 acc6 = _mm_mul_ps(_mm_load_ps(pa+20),_mm_load_ps(pb+20));
    __m128 acc7 = _mm_mul_ps(_mm_load_ps(pa+24),_mm_load_ps(pb+24));
    __m128 acc8 = _mm_mul_ps(_mm_load_ps(pa+28),_mm_load_ps(pb+28));
    pa += 32;
    pb += 32;
    while(--q!=0) { 
        __m128 i1 = _mm_mul_ps(_mm_load_ps(pa),_mm_load_ps(pb));
        __m128 i2 = _mm_mul_ps(_mm_load_ps(pa+4),_mm_load_ps(pb+4));
        __m128 i3 = _mm_mul_ps(_mm_load_ps(pa+8),_mm_load_ps(pb+8));
        __m128 i4 = _mm_mul_ps(_mm_load_ps(pa+12),_mm_load_ps(pb+12));
        __m128 i5 = _mm_mul_ps(_mm_load_ps(pa+16),_mm_load_ps(pb+16));
        __m128 i6 = _mm_mul_ps(_mm_load_ps(pa+20),_mm_load_ps(pb+20));
        __m128 i7 = _mm_mul_ps(_mm_load_ps(pa+24),_mm_load_ps(pb+24));
        __m128 i8 = _mm_mul_ps(_mm_load_ps(pa+28),_mm_load_ps(pb+28));
        pa += 32;
        pb += 32;
        acc1 = _mm_add_ps(acc1,i1);
        acc2 = _mm_add_ps(acc2,i2);
        acc3 = _mm_add_ps(acc3,i3);
        acc4 = _mm_add_ps(acc4,i4);
        acc5 = _mm_add_ps(acc5,i5);
        acc6 = _mm_add_ps(acc6,i6);
        acc7 = _mm_add_ps(acc7,i7);
        acc8 = _mm_add_ps(acc8,i8);
    }
    acc1 = _mm_add_ps(acc1,acc2);    
    acc3 = _mm_add_ps(acc3,acc4);    
    acc5 = _mm_add_ps(acc5,acc6);    
    acc7 = _mm_add_ps(acc7,acc8);    
    acc1 = _mm_add_ps(acc1,acc3);    
    acc5 = _mm_add_ps(acc5,acc7);    
    acc1 = _mm_add_ps(acc1,acc5);  
    acc1 = _mm_add_ps(acc1,_mm_movehl_ps(acc1,acc1));
    acc1 = _mm_add_ps(acc1,_mm_shuffle_ps(acc1,acc1,0x55));
    _mm_store_ss(&w,acc1);
    return w;
}

#else

ALGEBRA_INLINE float	vector_ps_float_aligned_32 (const float* pa,const float* pb,const size_t q) {
    return vector_ps_float_basic(pa,pb,q*32);
}

#endif


float	vector_ps_float (const float* v1,const float* v2,size_t n) {
    if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
        size_t q = n / 32;
        size_t r = n % 32;
        float w = 0;
        if (q > 0) w = vector_ps_float_aligned_32 (v1,v2,q);
        return w + vector_ps_float_basic (v1+q*32,v2+q*32,r);
    }
    return vector_ps_float_basic(v1,v2,n);
}

#define matrix_float_gran 1024

double	matrix_ps_float (const float* v1,const float* v2,size_t n,size_t m)
{
	size_t i;
	size_t q = (n*m) / matrix_float_gran;
	size_t r = (n*m) % matrix_float_gran;
    double s = 0;
    if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
		for (i=0;i<q;i++) {
			s += vector_ps_float_aligned_32 (v1+i*matrix_float_gran,v2+i*matrix_float_gran,matrix_float_gran/32);
		}
    }
	else {
		for (i=0;i<q;i++) {
		    s += vector_ps_float(v1+i*matrix_float_gran,v2+i*matrix_float_gran,matrix_float_gran);
		}
	}
	return s + vector_ps_float(v1+i*matrix_float_gran,v2+i*matrix_float_gran,r);
}

float	vector_sum_float_basic (const float* v,size_t n) 
{
    float w = 0;
    size_t i;
    for (i=0;i<n;i++)
        w += v[i];
    return w;
}

#if defined(ALGEBRA_SSE2)
 
ALGEBRA_INLINE float	vector_sum_float_aligned_32 (const float* pa,size_t q)
{
    float w;
    __m128 acc1 = _mm_load_ps(pa);
    __m128 acc2 = _mm_load_ps(pa+4);
    __m128 acc3 = _mm_load_ps(pa+8);
    __m128 acc4 = _mm_load_ps(pa+12); 
    __m128 acc5 = _mm_load_ps(pa+16); 
    __m128 acc6 = _mm_load_ps(pa+20);
    __m128 acc7 = _mm_load_ps(pa+24);
    __m128 acc8 = _mm_load_ps(pa+28); 
    pa += 32;
    while(--q!=0) {
        __m128 i1 = _mm_load_ps(pa);
        __m128 i2 = _mm_load_ps(pa+4);
        __m128 i3 = _mm_load_ps(pa+8);
        __m128 i4 = _mm_load_ps(pa+12); 
        __m128 i5 = _mm_load_ps(pa+16); 
        __m128 i6 = _mm_load_ps(pa+20);
        __m128 i7 = _mm_load_ps(pa+24);
        __m128 i8 = _mm_load_ps(pa+28); 
        pa += 32;
        acc1 = _mm_add_ps(acc1,i1);
        acc2 = _mm_add_ps(acc2,i2);
        acc3 = _mm_add_ps(acc3,i3);
        acc4 = _mm_add_ps(acc4,i4);
        acc5 = _mm_add_ps(acc5,i5);
        acc6 = _mm_add_ps(acc6,i6);
        acc7 = _mm_add_ps(acc7,i7);
        acc8 = _mm_add_ps(acc8,i8);
    }
    acc1 = _mm_add_ps(acc1,acc2);    
    acc3 = _mm_add_ps(acc3,acc4);    
    acc5 = _mm_add_ps(acc5,acc6);    
    acc7 = _mm_add_ps(acc7,acc8);    
    acc1 = _mm_add_ps(acc1,acc3);    
    acc5 = _mm_add_ps(acc5,acc7);    
    acc1 = _mm_add_ps(acc1,acc5);  
    acc1 = _mm_add_ps(acc1,_mm_movehl_ps(acc1,acc1));
    acc1 = _mm_add_ps(acc1,_mm_shuffle_ps(acc1,acc1,0x55));
    _mm_store_ss(&w,acc1);
    return w;
}

#else

ALGEBRA_INLINE float	vector_sum_float_aligned_32 (const float* pa,const size_t q) {
    return vector_sum_float_basic(pa,pb,q*32);
}

#endif


float	vector_sum_float (const float* v,size_t n) {
    if (ALGEBRA_IS_ALIGNED(v) ) {
        size_t q = n / 32;
        size_t r = n % 32;
        float w = 0;
        if (q > 0) w = vector_sum_float_aligned_32 (v,q);
        return w + vector_sum_float_basic (v+q*32,r);
    }
    return vector_sum_float_basic(v,n);
}

void	vector_cn2_float (float* z,float* x,float* y,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		z[i] = sqrt(x[i]*x[i]+y[i]*y[i]);
}

void	vector_qn2_float (float* z,float* x,float* yi,float* yj,float* yk,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		z[i] = sqrt(x[i]*x[i]+yi[i]*yi[i]+yj[i]*yj[i]+yk[i]*yk[i]);
}

float	vector_l1_float (float* v1,float* v2,size_t n)
{
	size_t i;
	float s,x;
	s = 0;
	for (i=0;i<n;i++)
	{
		x = v1[i] - v2[i];
		if (x > 0) s += x;
		else       s -= x;
	}
	return s;
}

float	vector_l2_float (const float* v1,const float* v2,size_t n)
{
    return sqrt(vector_l2p2_float(v1,v2,n));
}

#if defined(ALGEBRA_AVX)
float	vector_l2p2_float (const float* pa,const float* pb,size_t n) {

    
    size_t k;
    size_t q = n / 8;
    size_t r = n % 8;
    float w;
    if (q > 0) {
	    __m256 acc1 = _mm256_setzero_ps();
	    if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
	        for (k=0;k<q;k++) {
		        /* Charge 8 valeurs dans chaque tableau */
		        __m256 i1 = _mm256_load_ps(pa);
		        __m256 j1 = _mm256_load_ps(pb);
		        /* Différence */
		        __m256 d1 = _mm256_sub_ps(i1,j1);
		        /* Carré */
		        __m256 s1 = _mm256_mul_ps(d1,d1);
		        pa += 8;
		        pb += 8;
		        /* Accumule */
                acc1 = _mm256_add_ps(acc1,s1);
	        }
	    }
	    else {
	        for (k=0;k<q;k++) {
		        /* Charge 8 valeurs dans chaque tableau */
		        __m256 i1 = _mm256_loadu_ps(pa);
		        __m256 j1 = _mm256_loadu_ps(pb);
		        /* Différence */
		        __m256 d1 = _mm256_sub_ps(i1,j1);
		        /* Carré */
		        __m256 s1 = _mm256_mul_ps(d1,d1);
		        pa += 8;
		        pb += 8;
		        /* Accumule */
                acc1 = _mm256_add_ps(acc1,s1);
	        }
	    }

	    /* Somme finale */
        __m256 accp = _mm256_permute2f128_ps(acc1, acc1, 1);
        acc1 = _mm256_add_ps(acc1, accp);
        acc1 = _mm256_hadd_ps(acc1, acc1);
        acc1 = _mm256_hadd_ps(acc1, acc1);
        _mm_store_ss(&w,  _mm256_extractf128_ps(acc1,0));
    }
    else {
    	w = 0;
    }
    for (k=0;k<r;k++) {
	    float x = *pa++ - *pb++;
	    w += x*x;
    }
    return w;


}
#elif defined(ALGEBRA_SSSE3)

float	vector_l2p2_float (const float* pa,const float* pb,size_t n)
{
    size_t k;
    size_t q = n / 4;
    size_t r = n % 4;
    float w;
    if (q > 0) {
	__m128 acc1 = _mm_setzero_ps();
	if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
	    for (k=0;k<q;k++) {
		/* Charge 4 valeurs dans chaque tableau */
		__m128 i1 = _mm_load_ps(pa);
		__m128 j1 = _mm_load_ps(pb);
		/* Différence */
		__m128 d1 = _mm_sub_ps(i1,j1);
		/* Carré */
		__m128 s1 = _mm_mul_ps(d1,d1);
		pa += 4;
		pb += 4;
		/* Accumule */
                acc1 = _mm_add_ps(acc1,s1);
	    }
	}
	else {
	    for (k=0;k<q;k++) {
		/* Charge 4 valeurs dans chaque tableau */
		__m128 i1 = _mm_loadu_ps(pa);
		__m128 j1 = _mm_loadu_ps(pb);
		/* Différence */
		__m128 d1 = _mm_sub_ps(i1,j1);
		/* Carré */
		__m128 s1 = _mm_mul_ps(d1,d1);
		pa += 4;
		pb += 4;
		/* Accumule */
                acc1 = _mm_add_ps(acc1,s1);
	    }
	}
	/* Somme finale */
        acc1 = _mm_hadd_ps(acc1,acc1);
	_mm_store_ss(&w,_mm_hadd_ps(acc1,acc1));
    }
    else {
	w = 0;
    }
    for (k=0;k<r;k++) {
	float x = *pa++ - *pb++;
	w += x*x;
    }
    return w;
}

#else

float	vector_l2p2_float (const float* v1,const float* v2,size_t n)
{
	return vector_l2p2_float_basic(v1, v2, n);
}

#endif


float	vector_l2p2_float_basic (const float* v1,const float* v2,size_t n)
{
	size_t i;
	float s,x;
	s = 0;
	for (i=0;i<n;i++)
	{
		x = v1[i] - v2[i];
		s += x*x;
	}
	return s;
}


float	vector_chi1_float (float* v1,float* v2,size_t n)
{
	size_t i;
	float s,x,x1,y,minValue=1E-5;
	s = 0;
	for (i=0;i<n;i++)
	{
		x = v1[i] - v2[i];
		y = v1[i] + v2[i];
		x1 = -x;
		if (x1 > x) x = x1;
		if (y < minValue) y = minValue;
		s += x / y;
	}
	return s;
}

float	vector_chi2_float (const float* v1,const float* v2,size_t n)
{
	size_t i;
	float s,x,y,minValue=1E-5;
	s = 0;
	for (i=0;i<n;i++)
	{
		x = v1[i] - v2[i];
		y = v1[i] + v2[i];
		if (y < minValue) y = minValue;
		s += x * x / y;
	}
	return sqrt(s);
}

#ifdef ALGEBRA_SSSE3

float	vector_cos_float (const float* pa,const float* pb,size_t n)
{
    size_t k;
    size_t q = n / 4;
    size_t r = n % 4;
    double ps,norm,normA,normB;
    if (q > 0) {
	__m128d acc1 = _mm_setzero_pd();
	__m128d accA1 = _mm_setzero_pd();
	__m128d accB1 = _mm_setzero_pd();
	if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
	    for (k=0;k<q;k++) {
		/* Charge 4 valeurs dans chaque tableau */
		__m128 i1 = _mm_load_ps(pa);
		__m128 j1 = _mm_load_ps(pb);
		/* Multiplie */
		__m128 s1 = _mm_mul_ps(i1,j1); /* a*b */
                __m128 sA1 = _mm_mul_ps(i1,i1); /* a*a */
                __m128 sB1 = _mm_mul_ps(j1,j1); /* b*b */
		pa += 4;
		pb += 4;
		/* Accumule */
		acc1 = _mm_add_pd(acc1,_mm_cvtps_pd(_mm_hadd_ps(s1,s1)));
                accA1 = _mm_add_pd(accA1,_mm_cvtps_pd(_mm_hadd_ps(sA1,sA1)));
                accB1 = _mm_add_pd(accB1,_mm_cvtps_pd(_mm_hadd_ps(sB1,sB1)));
	    }
	}
	else {
	    for (k=0;k<q;k++) {
		/* Charge 4 valeurs dans chaque tableau */
		__m128 i1 = _mm_loadu_ps(pa);
		__m128 j1 = _mm_loadu_ps(pb);
		/* Multiplie */
		__m128 s1 = _mm_mul_ps(i1,j1); /* a*b */
                __m128 sA1 = _mm_mul_ps(i1,i1); /* a*a */
                __m128 sB1 = _mm_mul_ps(j1,j1); /* b*b */
		pa += 4;
		pb += 4;
		/* Accumule */
		acc1 = _mm_add_pd(acc1,_mm_cvtps_pd(_mm_hadd_ps(s1,s1)));
                accA1 = _mm_add_pd(accA1,_mm_cvtps_pd(_mm_hadd_ps(sA1,sA1)));
                accB1 = _mm_add_pd(accB1,_mm_cvtps_pd(_mm_hadd_ps(sB1,sB1)));
	    }
	}
	/* Somme finale */
	_mm_store_sd(&ps,_mm_hadd_pd(acc1,acc1));
	_mm_store_sd(&normA,_mm_hadd_pd(accA1,accA1));
	_mm_store_sd(&normB,_mm_hadd_pd(accB1,accB1));
    }
    else {
        ps = 0;
        normA = 0;
        normB = 0;
    }
    for (k=0;k<r;k++) {
	float a = *pa++;
	float b = *pb++;
	ps += a*b;
	normA += a*a;
	normB += b*b;
    }
    norm = sqrt(normA*normB);
    if (norm < 1E-5f)
	return 0;
    return ps / norm;
}

#else

float	vector_cos_float(const float* A,const float* B,size_t n)
{
    size_t i;
    double ps = 0;
    double normA = 0;
    double normB = 0;
    double norm;
    for (i=0;i<n;i++) {
	float a = A[i];
	float b = B[i];
	ps += a*b;
	normA += a*a;
	normB += b*b;
    }
    norm = sqrt(normA*normB);
    if (norm < 1E-5f)
	return 0;
    return ps / norm;
}

#endif

float	vector_linf_float (float* A,float* B,size_t n)
{
	size_t i;
	float s,x;
	s = A[0] - B[0];
	if (s < 0) s = -s;
	for (i=1;i<n;i++)
	{
		x = A[i] - B[i];
		if (x < 0) x = -x;
		if (x > s) x = s;
	}
	return s;
}

float	vector_std_float (const float* v1,const float* v2,size_t n)
{
    size_t i;
    float mean = 0, var = 0;
    for (i=0;i<n;i++) {
        float x = (v1[i]-v2[i]);
        mean += x;
        var += x*x;
    }
    mean /= n;
    return sqrt(var/n-mean*mean);
}

void	vector_cpym_float (float* v1,float lambda,const float* v2,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v1[i] = lambda*v2[i];
}

#ifdef ALGEBRA_SSE2

void	vector_add_float (float* v1,const float* v2,size_t n)
{
	size_t k;
	
	size_t q = n / 16;
	size_t r = n % 16;
	if(q > 0) {
		if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
			for (k=0;k<q;k++) {
				/* Charge 4 valeurs de chaque tableau */
				__m128 i1 = _mm_load_ps(v1);
				__m128 j1 = _mm_load_ps(v2);
				__m128 i2 = _mm_load_ps(v1+4);
				__m128 j2 = _mm_load_ps(v2+4);
				__m128 i3 = _mm_load_ps(v1+8);
				__m128 j3 = _mm_load_ps(v2+8);
				__m128 i4 = _mm_load_ps(v1+12);
				__m128 j4 = _mm_load_ps(v2+12);
				/* Multiplie */
				i1 = _mm_add_ps(i1,j1);
				i2 = _mm_add_ps(i2,j2);
				i3 = _mm_add_ps(i3,j3);
				i4 = _mm_add_ps(i4,j4);
				/* Sauvegarde */
				_mm_store_ps(v1, i1);
				_mm_store_ps(v1+4, i2);
				_mm_store_ps(v1+8, i3);
				_mm_store_ps(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
		else {		
			for (k=0;k<q;k++) {
				/* Charge 4 valeurs de chaque tableau */
				__m128 i1 = _mm_loadu_ps(v1);
				__m128 j1 = _mm_loadu_ps(v2);
				__m128 i2 = _mm_loadu_ps(v1+4);
				__m128 j2 = _mm_loadu_ps(v2+4);
				__m128 i3 = _mm_loadu_ps(v1+8);
				__m128 j3 = _mm_loadu_ps(v2+8);
				__m128 i4 = _mm_loadu_ps(v1+12);
				__m128 j4 = _mm_loadu_ps(v2+12);
				/* Multiplie */
				i1 = _mm_add_ps(i1,j1);
				i2 = _mm_add_ps(i2,j2);
				i3 = _mm_add_ps(i3,j3);
				i4 = _mm_add_ps(i4,j4);
				/* Sauvegarde */
				_mm_storeu_ps(v1, i1);
				_mm_storeu_ps(v1+4, i2);
				_mm_storeu_ps(v1+8, i3);
				_mm_storeu_ps(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
	}
	
	for(k = 0 ; k<r ; k++)
		v1[k] += v2[k];
}

#else

void	vector_add_float (float* v1,const float* v2,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v1[i] += v2[i];
}

#endif


#if defined(ALGEBRA_AVX)

ALGEBRA_INLINE void		vector_addm_float_aligned_32 (float* v1,float lambda,const float* v2,size_t n)
{
	size_t k;
	
	__m256 l1 = _mm256_broadcast_ss(&lambda);

	size_t q = n / 32;
	size_t r = n % 32;
	if(q > 0) {
		if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
			for (k=0;k<q;k++) {
				/* Charge 8 valeurs de chaque tableau */
				__m256 i1 = _mm256_load_ps(v1);
				__m256 j1 = _mm256_load_ps(v2);
					   j1 = _mm256_mul_ps(j1, l1);
				__m256 i2 = _mm256_load_ps(v1+8);
				__m256 j2 = _mm256_load_ps(v2+8);
					   j2 = _mm256_mul_ps(j2, l1);
				__m256 i3 = _mm256_load_ps(v1+16);
				__m256 j3 = _mm256_load_ps(v2+16);
					   j3 = _mm256_mul_ps(j3, l1);
				__m256 i4 = _mm256_load_ps(v1+24);
				__m256 j4 = _mm256_load_ps(v2+24);
					   j4 = _mm256_mul_ps(j4, l1);
				/* Soustrait */
				i1 = _mm256_add_ps(i1,j1);
				i2 = _mm256_add_ps(i2,j2);
				i3 = _mm256_add_ps(i3,j3);
				i4 = _mm256_add_ps(i4,j4);
				/* Sauvegarde */
				_mm256_store_ps(v1, i1);
				_mm256_store_ps(v1+8, i2);
				_mm256_store_ps(v1+16, i3);
				_mm256_store_ps(v1+24, i4);
				v1 += 32;
				v2 += 32;
			}
		}
		else {		
			for (k=0;k<q;k++) {
				/* Charge 8 valeurs de chaque tableau */
				__m256 i1 = _mm256_loadu_ps(v1);
				__m256 j1 = _mm256_loadu_ps(v2);
					   j1 = _mm256_mul_ps(j1, l1);
				__m256 i2 = _mm256_loadu_ps(v1+8);
				__m256 j2 = _mm256_loadu_ps(v2+8);
					   j2 = _mm256_mul_ps(j2, l1);
				__m256 i3 = _mm256_loadu_ps(v1+16);
				__m256 j3 = _mm256_loadu_ps(v2+16);
					   j3 = _mm256_mul_ps(j3, l1);
				__m256 i4 = _mm256_loadu_ps(v1+24);
				__m256 j4 = _mm256_loadu_ps(v2+24);
					   j4 = _mm256_mul_ps(j4, l1);
				/* Soustrait */
				i1 = _mm256_add_ps(i1,j1);
				i2 = _mm256_add_ps(i2,j2);
				i3 = _mm256_add_ps(i3,j3);
				i4 = _mm256_add_ps(i4,j4);
				/* Sauvegarde */
				_mm256_storeu_ps(v1, i1);
				_mm256_storeu_ps(v1+8, i2);
				_mm256_storeu_ps(v1+16, i3);
				_mm256_storeu_ps(v1+24, i4);
				v1 += 32;
				v2 += 32;
			}
		}
	}
	
	for(k = 0 ; k<r ; k++)
		v1[k] += lambda*v2[k];
}


#elif defined(ALGEBRA_SSSE3)

ALGEBRA_INLINE void		vector_addm_float_aligned_32 (float* v1,float lambda,const float* v2,size_t n)
{
	size_t k;
	
	__m128 l1 = _mm_load1_ps(&lambda);

	size_t q = n / 16;
	size_t r = n % 16;
	if(q > 0) {
		if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
			for (k=0;k<q;k++) {
				/* Charge 8 valeurs de chaque tableau */
				__m128 i1 = _mm_load_ps(v1);
				__m128 j1 = _mm_load_ps(v2);
					   j1 = _mm_mul_ps(j1, l1);
				__m128 i2 = _mm_load_ps(v1+4);
				__m128 j2 = _mm_load_ps(v2+4);
					   j2 = _mm_mul_ps(j2, l1);
				__m128 i3 = _mm_load_ps(v1+8);
				__m128 j3 = _mm_load_ps(v2+8);
					   j3 = _mm_mul_ps(j3, l1);
				__m128 i4 = _mm_load_ps(v1+12);
				__m128 j4 = _mm_load_ps(v2+12);
					   j4 = _mm_mul_ps(j4, l1);
				/* Soustrait */
				i1 = _mm_add_ps(i1,j1);
				i2 = _mm_add_ps(i2,j2);
				i3 = _mm_add_ps(i3,j3);
				i4 = _mm_add_ps(i4,j4);
				/* Sauvegarde */
				_mm_store_ps(v1, i1);
				_mm_store_ps(v1+4, i2);
				_mm_store_ps(v1+8, i3);
				_mm_store_ps(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
		else {		
			for (k=0;k<q;k++) {
				/* Charge 8 valeurs de chaque tableau */
				__m128 i1 = _mm_loadu_ps(v1);
				__m128 j1 = _mm_loadu_ps(v2);
					   j1 = _mm_mul_ps(j1, l1);
				__m128 i2 = _mm_loadu_ps(v1+4);
				__m128 j2 = _mm_loadu_ps(v2+4);
					   j2 = _mm_mul_ps(j2, l1);
				__m128 i3 = _mm_loadu_ps(v1+8);
				__m128 j3 = _mm_loadu_ps(v2+8);
					   j3 = _mm_mul_ps(j3, l1);
				__m128 i4 = _mm_loadu_ps(v1+12);
				__m128 j4 = _mm_loadu_ps(v2+12);
					   j4 = _mm_mul_ps(j4, l1);
				/* Soustrait */
				i1 = _mm_add_ps(i1,j1);
				i2 = _mm_add_ps(i2,j2);
				i3 = _mm_add_ps(i3,j3);
				i4 = _mm_add_ps(i4,j4);
				/* Sauvegarde */
				_mm_storeu_ps(v1, i1);
				_mm_storeu_ps(v1+4, i2);
				_mm_storeu_ps(v1+8, i3);
				_mm_storeu_ps(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
	}
	
	for(k = 0 ; k<r ; k++)
		v1[k] += lambda*v2[k];

}

#else


ALGEBRA_INLINE void	vector_addm_float_aligned_32 (float* v1,float lambda,const float* v2,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v1[i] += lambda*v2[i];
}

#endif

ALGEBRA_INLINE void	vector_addm_float_basic (float* v1,float lambda,const float* v2,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v1[i] += lambda*v2[i];
}


void	vector_addm_float (float* v1,float lambda,const float* v2,size_t n)
{
	vector_addm_float_aligned_32(v1, lambda, v2, n);
}

void	vector_sub_float_basic (float* v1,const float* v2,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v1[i] -= v2[i];
}


#if defined(ALGEBRA_AVX)

void	vector_sub_float (float* v1,const float* v2,size_t n)
{
	size_t k;

	size_t q = n / 32;
	size_t r = n % 32;
	if(q > 0) {
		if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
			for (k=0;k<q;k++) {
				/* Charge 8 valeurs de chaque tableau */
				__m256 i1 = _mm256_load_ps(v1);
				__m256 j1 = _mm256_load_ps(v2);
				__m256 i2 = _mm256_load_ps(v1+8);
				__m256 j2 = _mm256_load_ps(v2+8);
				__m256 i3 = _mm256_load_ps(v1+16);
				__m256 j3 = _mm256_load_ps(v2+16);
				__m256 i4 = _mm256_load_ps(v1+24);
				__m256 j4 = _mm256_load_ps(v2+24);
				/* Soustrait */
				i1 = _mm256_sub_ps(i1,j1);
				i2 = _mm256_sub_ps(i2,j2);
				i3 = _mm256_sub_ps(i3,j3);
				i4 = _mm256_sub_ps(i4,j4);
				/* Sauvegarde */
				_mm256_store_ps(v1, i1);
				_mm256_store_ps(v1+8, i2);
				_mm256_store_ps(v1+16, i3);
				_mm256_store_ps(v1+24, i4);
				v1 += 32;
				v2 += 32;
			}
		}
		else {		
			for (k=0;k<q;k++) {
				/* Charge 8 valeurs de chaque tableau */
				__m256 i1 = _mm256_loadu_ps(v1);
				__m256 j1 = _mm256_loadu_ps(v2);
				__m256 i2 = _mm256_loadu_ps(v1+8);
				__m256 j2 = _mm256_loadu_ps(v2+8);
				__m256 i3 = _mm256_loadu_ps(v1+16);
				__m256 j3 = _mm256_loadu_ps(v2+16);
				__m256 i4 = _mm256_loadu_ps(v1+24);
				__m256 j4 = _mm256_loadu_ps(v2+24);
				/* Soustrait */
				i1 = _mm256_sub_ps(i1,j1);
				i2 = _mm256_sub_ps(i2,j2);
				i3 = _mm256_sub_ps(i3,j3);
				i4 = _mm256_sub_ps(i4,j4);
				/* Sauvegarde */
				_mm256_storeu_ps(v1, i1);
				_mm256_storeu_ps(v1+8, i2);
				_mm256_storeu_ps(v1+16, i3);
				_mm256_storeu_ps(v1+24, i4);
				v1 += 32;
				v2 += 32;
			}
		}
	}
	
	for(k = 0 ; k<r ; k++)
		v1[k] -= v2[k];
}


#elif defined(ALGEBRA_SSE2)

void	vector_sub_float (float* v1,const float* v2,size_t n)
{
	size_t k;
	
	size_t q = n / 16;
	size_t r = n % 16;
	if(q > 0) {
		if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
			for (k=0;k<q;k++) {
				/* Charge 4 valeurs de chaque tableau */
				__m128 i1 = _mm_load_ps(v1);
				__m128 j1 = _mm_load_ps(v2);
				__m128 i2 = _mm_load_ps(v1+4);
				__m128 j2 = _mm_load_ps(v2+4);
				__m128 i3 = _mm_load_ps(v1+8);
				__m128 j3 = _mm_load_ps(v2+8);
				__m128 i4 = _mm_load_ps(v1+12);
				__m128 j4 = _mm_load_ps(v2+12);
				/* Multiplie */
				i1 = _mm_sub_ps(i1,j1);
				i2 = _mm_sub_ps(i2,j2);
				i3 = _mm_sub_ps(i3,j3);
				i4 = _mm_sub_ps(i4,j4);
				/* Sauvegarde */
				_mm_store_ps(v1, i1);
				_mm_store_ps(v1+4, i2);
				_mm_store_ps(v1+8, i3);
				_mm_store_ps(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
		else {		
			for (k=0;k<q;k++) {
				/* Charge 4 valeurs de chaque tableau */
				__m128 i1 = _mm_loadu_ps(v1);
				__m128 j1 = _mm_loadu_ps(v2);
				__m128 i2 = _mm_loadu_ps(v1+4);
				__m128 j2 = _mm_loadu_ps(v2+4);
				__m128 i3 = _mm_loadu_ps(v1+8);
				__m128 j3 = _mm_loadu_ps(v2+8);
				__m128 i4 = _mm_loadu_ps(v1+12);
				__m128 j4 = _mm_loadu_ps(v2+12);
				/* Multiplie */
				i1 = _mm_sub_ps(i1,j1);
				i2 = _mm_sub_ps(i2,j2);
				i3 = _mm_sub_ps(i3,j3);
				i4 = _mm_sub_ps(i4,j4);
				/* Sauvegarde */
				_mm_storeu_ps(v1, i1);
				_mm_storeu_ps(v1+4, i2);
				_mm_storeu_ps(v1+8, i3);
				_mm_storeu_ps(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
	}
	
	for(k = 0 ; k<r ; k++)
		v1[k] -= v2[k];
}

#else

void	vector_sub_float (float* v1,const float* v2,size_t n)
{
	vector_sub_float_basic(v1, v2, n);
}

#endif

#ifdef ALGEBRA_SSE2

void	vector_mul_float (float* v1,float* v2,size_t n)
{
	size_t k;
	
	size_t q = n / 16;
	size_t r = n % 16;
	if(q > 0) {
		if (ALGEBRA_IS_ALIGNED(v1) && ALGEBRA_IS_ALIGNED(v2)) {
			for (k=0;k<q;k++) {
				/* Charge 4 valeurs de chaque tableau */
				__m128 i1 = _mm_load_ps(v1);
				__m128 j1 = _mm_load_ps(v2);
				__m128 i2 = _mm_load_ps(v1+4);
				__m128 j2 = _mm_load_ps(v2+4);
				__m128 i3 = _mm_load_ps(v1+8);
				__m128 j3 = _mm_load_ps(v2+8);
				__m128 i4 = _mm_load_ps(v1+12);
				__m128 j4 = _mm_load_ps(v2+12);
				/* Multiplie */
				i1 = _mm_mul_ps(i1,j1);
				i2 = _mm_mul_ps(i2,j2);
				i3 = _mm_mul_ps(i3,j3);
				i4 = _mm_mul_ps(i4,j4);
				/* Sauvegarde */
				_mm_store_ps(v1, i1);
				_mm_store_ps(v1+4, i2);
				_mm_store_ps(v1+8, i3);
				_mm_store_ps(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
		else {		
			for (k=0;k<q;k++) {
				/* Charge 4 valeurs de chaque tableau */
				__m128 i1 = _mm_loadu_ps(v1);
				__m128 j1 = _mm_loadu_ps(v2);
				__m128 i2 = _mm_loadu_ps(v1+4);
				__m128 j2 = _mm_loadu_ps(v2+4);
				__m128 i3 = _mm_loadu_ps(v1+8);
				__m128 j3 = _mm_loadu_ps(v2+8);
				__m128 i4 = _mm_loadu_ps(v1+12);
				__m128 j4 = _mm_loadu_ps(v2+12);
				/* Multiplie */
				i1 = _mm_mul_ps(i1,j1);
				i2 = _mm_mul_ps(i2,j2);
				i3 = _mm_mul_ps(i3,j3);
				i4 = _mm_mul_ps(i4,j4);
				/* Sauvegarde */
				_mm_storeu_ps(v1, i1);
				_mm_storeu_ps(v1+4, i2);
				_mm_storeu_ps(v1+8, i3);
				_mm_storeu_ps(v1+12, i4);
				v1 += 16;
				v2 += 16;
			}
		}
	}
	
	for(k = 0 ; k<r ; k++)
		v1[k] *= v2[k];
}

#else

void	vector_mul_float (float* v1,float* v2,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v1[i] *= v2[i];
}

#endif

void	vector_div_float (float* v1,float* v2,size_t n)
{
	size_t i;
	for (i=0;i<n;i++)
		v1[i] /= v2[i];
}


size_t  vector_argmin_l2_float (const float* x,const float* C,size_t n,size_t q,const float* normC)
{
    if (ALGEBRA_IS_ALIGNED(x) && ALGEBRA_IS_ALIGNED(C) && ALGEBRA_IS_ALIGNED(C+n) && (n%32) == 0) {
        size_t p = n / 32;
        size_t c,cMin = 0;
        float dMin = 1E37;
        for (c=0;c<q;c++) {
            float d = normC[c] - 2*vector_ps_float_aligned_32(x,C+c*n,p);
            if (d < dMin) {
                dMin = d;
                cMin = c;
            }
        }
        return cMin;
    }
    else {
        /*printf ("vector_argmin_l2_float: unaligned data");*/
        return ~0;
    }    
}

int     vector_isfinite_float (const float* v,size_t n)
{
    size_t i;
    for (i=0;i<n;i++) {
        if (!isfinite(v[i]))
            return 0;
    }
    return 1;
}
