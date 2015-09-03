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
 * \file core.h
 * \author Philippe H. Gosselin
 * \version 1.0
 */

#ifndef __algebra_core_h__
#define __algebra_core_h__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <malloc.h>

#ifdef __SSE2__
#define ALGEBRA_SSE2
#include <emmintrin.h>
void _mm_print_ps(__m128 x);

#ifdef __SSSE3__
#define ALGEBRA_SSSE3
#include <tmmintrin.h>

#ifdef __SSE4_1__
#define ALGEBRA_SSE4
#include <smmintrin.h>

#ifdef __AVX__
#define ALGEBRA_AVX
#include <immintrin.h>
void _mm256_print_ps(__m256 x);

#endif /* AVX */

#endif /* SSE4_1 */

#endif /* SSSE3 */

#ifdef __AVX__
#define ALGEBRA_FLOAT_GRANULARITY   8
#define ALGEBRA_ALIGN_SIZE          0x20
#else /* Not AVX */
#define ALGEBRA_FLOAT_GRANULARITY   4
#define ALGEBRA_ALIGN_SIZE          0x10
#endif

#define ALGEBRA_IS_ALIGNED(x) ((((size_t)x)%ALGEBRA_ALIGN_SIZE) == 0)

#define algebra_alloc_uchar(n)  ((unsigned char*)_mm_malloc((n)*sizeof(unsigned char),ALGEBRA_ALIGN_SIZE))
#define algebra_alloc_short(n)  ((short*)_mm_malloc((n)*sizeof(short),ALGEBRA_ALIGN_SIZE))
#define algebra_alloc_float(n)  ((float*)_mm_malloc((n)*sizeof(float),ALGEBRA_ALIGN_SIZE))
#define algebra_alloc_double(n) ((double*)_mm_malloc((n)*sizeof(double),ALGEBRA_ALIGN_SIZE))

#define algebra_free(p)   _mm_free(p)

#ifdef __cplusplus
template<class T>
T* algebra_alloc (size_t n) {
	return (T*)_mm_malloc(n*sizeof(T),ALGEBRA_ALIGN_SIZE);
}
#endif

#else /* Not SSE2 */

#define ALGEBRA_IS_ALIGNED(x) (1)

#define algebra_alloc_uchar(n)  ((unsigned char*)malloc((n)*sizeof(unsigned char)))
#define algebra_alloc_short(n)  ((short*)malloc(n*sizeof(short)))
#define algebra_alloc_float(n)  ((float*)malloc(n*sizeof(float)))
#define algebra_alloc_double(n) ((double*)malloc(n*sizeof(double)))

#define algebra_free(p)   free(p)

#ifdef __cplusplus
template<class T>
T* algebra_alloc(size_t n) {
	return (T*) malloc(n * sizeof(T));
}
#endif

#endif /* SSE2 */

#if defined(__STRICT_ANSI__)  
#define ALGEBRA_INLINE 
#elif defined(__GNUC__)
#define ALGEBRA_INLINE inline
#else
#define ALGEBRA_INLINE 
#endif

#ifdef __cplusplus
extern "C" {
#endif

void algebra_show_arch();

#ifdef __cplusplus
}
#endif

#endif
