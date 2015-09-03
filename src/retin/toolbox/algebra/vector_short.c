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

#include "vector_short.h"

#include <stdio.h>


#ifdef ALGEBRA_SSE4

int	vector_ps_short (const short* pa,const short* pb,size_t n)
{
    size_t k;
    size_t q = n / 16;
    size_t r = n % 16;
    int w;
    if (q > 0) {
	__m128i acc1 = _mm_setzero_si128();
	__m128i acc2 = _mm_setzero_si128();
	if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
	    for (k=0;k<q;k++) {
		/* Charge 16 mots dans chaque tableau */
		__m128i a1 = _mm_load_si128((__m128i*)pa);
		__m128i b1 = _mm_load_si128((__m128i*)pb);
		__m128i a2 = _mm_load_si128((__m128i*)(pa+8));
		__m128i b2 = _mm_load_si128((__m128i*)(pb+8));
		/* Multiple, somme et converti en double word */
		__m128i s1 = _mm_madd_epi16(a1,b1);
		__m128i s2 = _mm_madd_epi16(a2,b2);
		pa += 16;
		pb += 16;
		/* Accumule */
		acc1 = _mm_add_epi32(acc1,s1);
		acc2 = _mm_add_epi32(acc2,s2);
	    }
	}
	else {
	    for (k=0;k<q;k++) {
		/* Charge 16 mots dans chaque tableau */
		__m128i a1 = _mm_loadu_si128((__m128i*)pa);
		__m128i b1 = _mm_loadu_si128((__m128i*)pb);
		__m128i a2 = _mm_loadu_si128((__m128i*)(pa+8));
		__m128i b2 = _mm_loadu_si128((__m128i*)(pb+8));
		/* Multiple, somme et converti en double word */
		__m128i s1 = _mm_madd_epi16(a1,b1);
		__m128i s2 = _mm_madd_epi16(a2,b2);
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

int	vector_ps_short (const short* v1,const short* v2,size_t n)
{
	size_t i;
	int s = 0;
	for (i=0;i<n;i++)
		s += v1[i]*v2[i];
	return s;
}

#endif

#ifdef ALGEBRA_SSE4

float	vector_cos_short (const short* pa,const short* pb,size_t n)
{
    size_t k;
    double norm;
    size_t q = n / 16;
    size_t r = n % 16;
    int ps,na,nb;
    if (q > 0) {
        __m128i acc;
	__m128i acc_ps1 = _mm_setzero_si128();
	__m128i acc_ps2 = _mm_setzero_si128();
	__m128i acc_na1 = _mm_setzero_si128();
	__m128i acc_na2 = _mm_setzero_si128();
	__m128i acc_nb1 = _mm_setzero_si128();
	__m128i acc_nb2 = _mm_setzero_si128();
	if (ALGEBRA_IS_ALIGNED(pa) && ALGEBRA_IS_ALIGNED(pb)) {
	    for (k=0;k<q;k++) {
		/* Charge 16 mots dans chaque tableau */
		__m128i a1 = _mm_load_si128((__m128i*)pa);
		__m128i b1 = _mm_load_si128((__m128i*)pb);
		__m128i a2 = _mm_load_si128((__m128i*)(pa+8));
		__m128i b2 = _mm_load_si128((__m128i*)(pb+8));
		/* Multiple, somme et converti en double word */
		__m128i ps1 = _mm_madd_epi16(a1,b1);
		__m128i ps2 = _mm_madd_epi16(a2,b2);
		__m128i na1 = _mm_madd_epi16(a1,a1);
		__m128i na2 = _mm_madd_epi16(a2,a2);
		__m128i nb1 = _mm_madd_epi16(b1,b1);
		__m128i nb2 = _mm_madd_epi16(b2,b2);
		pa += 16;
		pb += 16;
		/* Accumule */
		acc_ps1 = _mm_add_epi32(acc_ps1,ps1);
		acc_ps2 = _mm_add_epi32(acc_ps2,ps2);
		acc_na1 = _mm_add_epi32(acc_na1,na1);
		acc_na2 = _mm_add_epi32(acc_na2,na2);
		acc_nb1 = _mm_add_epi32(acc_nb1,nb1);
		acc_nb2 = _mm_add_epi32(acc_nb2,nb2);
	    }
	}
	else {
	    for (k=0;k<q;k++) {
	    }
	}
	/* Somme finale */
	acc = _mm_add_epi32(acc_ps1,acc_ps2);
	acc = _mm_hadd_epi32(acc,acc);
	acc = _mm_hadd_epi32(acc,acc);
	ps = _mm_extract_epi32(acc,0);

	acc = _mm_add_epi32(acc_na1,acc_na2);
	acc = _mm_hadd_epi32(acc,acc);
	acc = _mm_hadd_epi32(acc,acc);
	na = _mm_extract_epi32(acc,0);

	acc = _mm_add_epi32(acc_nb1,acc_nb2);
	acc = _mm_hadd_epi32(acc,acc);
	acc = _mm_hadd_epi32(acc,acc);
	nb = _mm_extract_epi32(acc,0);
    }
    else {
	ps = 0;
	na = 0;
	nb = 0;
    }
    for (k=0;k<r;k++) {
	int a = *pa++;
	int b = *pb++;
	ps += a*b;
	na += a*a;
	nb += b*b;
    }
    norm = sqrt( ((double)na) * ((double)nb) );
    if (norm < 1E-5f)
	return 0;
    return ps / norm;
}

#else

float	vector_cos_short(const short* A,const short* B,size_t n)
{
    size_t i;
    int ps = 0;
    int normA = 0;
    int normB = 0;
    double norm;
    for (i=0;i<n;i++) {
	int a = A[i];
	int b = B[i];
	ps += a*b;
	normA += a*a;
	normB += b*b;
    }
    norm = sqrt( ((double)normA) * ((double)normB) );
    if (norm < 1E-5f)
	return 0;
    return ps / norm;
}

#endif
