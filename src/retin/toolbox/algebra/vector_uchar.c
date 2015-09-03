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

#include "vector_uchar.h"

#include <stdio.h>
 
#ifdef ALGEBRA_SSE4

unsigned int	vector_ps_uchar (const unsigned char* pa,const unsigned char* pb,size_t n)
{
    size_t k;
    size_t q = n / 16;
    size_t r = n % 16;
    unsigned int w;
    if (q > 0) {
	__m128i zero = _mm_setzero_si128();
	__m128i acc1 = _mm_setzero_si128();
	__m128i acc2 = _mm_setzero_si128();
	if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
	    for (k=0;k<q;k++) {
		/* Charge 16 octets dans chaque tableau */
		__m128i a = _mm_load_si128((__m128i*)pa);
		__m128i b = _mm_load_si128((__m128i*)pb);
		/* Transforme en word */
		__m128i s1a = _mm_cvtepu8_epi16(a);
		__m128i s2a = _mm_unpackhi_epi8(a,zero);
		__m128i s1b = _mm_cvtepu8_epi16(b);
		__m128i s2b = _mm_unpackhi_epi8(b,zero);
		/* Multiple, somme et converti en double word */
		__m128i s1 = _mm_madd_epi16(s1a,s1b);
		__m128i s2 = _mm_madd_epi16(s2a,s2b);
		pa += 16;
		pb += 16;
		/* Accumule */
		acc1 = _mm_add_epi32(acc1,s1);
		acc2 = _mm_add_epi32(acc2,s2);
	    }
	}
	else {
	    for (k=0;k<q;k++) {
		/* Charge 16 octets dans chaque tableau */
		__m128i a = _mm_loadu_si128((__m128i*)pa);
		__m128i b = _mm_loadu_si128((__m128i*)pb);
		/* Transforme en word */
		__m128i s1a = _mm_unpacklo_epi8(a,zero);
		__m128i s2a = _mm_unpackhi_epi8(a,zero);
		__m128i s1b = _mm_unpacklo_epi8(b,zero);
		__m128i s2b = _mm_unpackhi_epi8(b,zero);
		/* Multiple, somme et converti en double word */
		__m128i s1 = _mm_madd_epi16(s1a,s1b);
		__m128i s2 = _mm_madd_epi16(s2a,s2b);
		pa += 16;
		pb += 16;
		/* Accumule */
		acc1 = _mm_add_epi32(acc1,s1);
		acc2 = _mm_add_epi32(acc2,s2);
	    }
	}
	/* Somme finale */
	acc1 = _mm_add_epi32(acc1,acc2);
	acc1 = _mm_hadd_epi32(acc1,acc1);
	acc1 = _mm_hadd_epi32(acc1,acc1);
	w = _mm_extract_epi32(acc1,0);
    }
    else {
	w = 0;
    }
    for (k=0;k<r;k++)
	w += (*pa++) * (*pb++);
    return w;
}

#else

unsigned int	vector_ps_uchar (const unsigned char* v1,const unsigned char* v2,size_t n)
{
	size_t i;
	unsigned int s = 0;
	for (i=0;i<n;i++)
		s += v1[i]*v2[i];
	return s;
}

#endif

#ifdef ALGEBRA_SSE4

unsigned int	vector_l1_uchar (const unsigned char* pa,const unsigned char* pb,size_t n)
{
    size_t k;
    size_t q = n / 16;
    size_t r = n % 16;
    unsigned int w;
    if (q > 0) {
	__m128i acc = _mm_setzero_si128();
	if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
	    for (k=0;k<q;k++) {
		/* Charge 16 octets dans chaque tableau */
		__m128i i1 = _mm_load_si128((__m128i*)pa);
		__m128i j1 = _mm_load_si128((__m128i*)pb);
		/* Somme des valeurs absolues de la différence entre les 16 bytes */
		__m128i s1 = _mm_sad_epu8(i1,j1);
		pa += 16;
		pb += 16;
		/* Accumule */
		acc = _mm_add_epi32(acc,s1);
	    }
	}
	else {
	    for (k=0;k<q;k++) {
		/* Charge 16 octets dans chaque tableau */
		__m128i i1 = _mm_loadu_si128((__m128i*)pa);
		__m128i j1 = _mm_loadu_si128((__m128i*)pb);
		/* Somme des valeurs absolues de la différence entre les 16 bytes */
		__m128i s1 = _mm_sad_epu8(i1,j1);
		pa += 16;
		pb += 16;
		/* Accumule */
		acc = _mm_add_epi32(acc,s1);
	    }
	}
	/* Somme finale */
	acc = _mm_hadd_epi32(acc,acc);
	acc = _mm_hadd_epi32(acc,acc);
 	w = _mm_cvtsi128_si32(acc);
    }
    else {
	w = 0;
    }
    for (k=0;k<r;k++) {
	int x1 = *pa++ - *pb++;
	int x2 = -x1;
	w += (x1>x2)?x1:x2;
    }
    return w;
}

#else

unsigned int	vector_l1_uchar (const unsigned char* pa,const unsigned char* pb,size_t n)
{
    size_t k;
    unsigned int w = 0;
    for (k=0;k<n;k++) {
	int x1 = *pa++ - *pb++;
	int x2 = -x1;
	w += (x1>x2)?x1:x2;
    }
    return w;
}

#endif

#ifdef ALGEBRA_SSE4

unsigned int	vector_l2p2_uchar (const unsigned char* pa,const unsigned char* pb,size_t n)
{
    size_t k;
    size_t q = n / 16;
    size_t r = n % 16;
    unsigned int x,w;
    if (q > 0) {
	__m128i zero = _mm_setzero_si128();
	__m128i acc1 = _mm_setzero_si128();
	__m128i acc2 = _mm_setzero_si128();
	if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
	    for (k=0;k<q;k++) {
		/* Charge 16 octets dans chaque tableau */
		__m128i a = _mm_load_si128((__m128i*)pa);
		__m128i b = _mm_load_si128((__m128i*)pb);
		/* Valeur absolue de la différence entre les 16 bytes */
		__m128i d = _mm_or_si128(_mm_subs_epu8(a,b),_mm_subs_epu8(b,a));
		/* Transforme en word */
		__m128i s1 = _mm_unpacklo_epi8(d,zero);
		__m128i s2 = _mm_unpackhi_epi8(d,zero);
		pa += 16;
		pb += 16;
		/* Monte au carré, somme et converti en double word */
		s1 = _mm_madd_epi16(s1,s1);
		s2 = _mm_madd_epi16(s2,s2);
		/* Accumule */
		acc1 = _mm_add_epi32(acc1,s1);
		acc2 = _mm_add_epi32(acc2,s2);
	    }
	}
	else {
	    for (k=0;k<q;k++) {
		/* Charge 16 octets dans chaque tableau */
		__m128i a = _mm_loadu_si128((__m128i*)pa);
		__m128i b = _mm_loadu_si128((__m128i*)pb);
		/* Valeur absolue de la différence entre les 16 bytes */
		__m128i d = _mm_or_si128(_mm_subs_epu8(a,b),_mm_subs_epu8(b,a));
		/* Transforme en word */
		__m128i s1 = _mm_unpacklo_epi8(d,zero);
		__m128i s2 = _mm_unpackhi_epi8(d,zero);
		pa += 16;
		pb += 16;
		/* Monte au carré, somme et converti en double word */
		s1 = _mm_madd_epi16(s1,s1);
		s2 = _mm_madd_epi16(s2,s2);
		/* Accumule */
		acc1 = _mm_add_epi32(acc1,s1);
		acc2 = _mm_add_epi32(acc2,s2);
	    }
	}
	/* Somme finale */
	acc1 = _mm_add_epi32(acc1,acc2);
	acc1 = _mm_hadd_epi32(acc1,acc1);
	acc1 = _mm_hadd_epi32(acc1,acc1);
	w = _mm_extract_epi32(acc1,0);
    }
    else {
	w = 0;
    }
    for (k=0;k<r;k++) {
	x = (*pa++) - (*pb++);
	w += x*x;
    }
    return w;
}

#else

unsigned int	vector_l2p2_uchar (const unsigned char* v1,const unsigned char* v2,size_t n)
{
	size_t i;
	unsigned int x,s = 0;
	for (i=0;i<n;i++)
	{
		x = v1[i] - v2[i];
		s += x*x;
	}
	return s;
}

#endif
