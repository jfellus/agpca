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
 * \file core.c
 * \author Philippe H. Gosselin
 * \version 1.0
 */

#include "core.h"

#include <stdio.h>
 
#ifdef ALGEBRA_SSE2

void _mm_print_ps(__m128 x) {
    size_t i;
    float tab[4];
    _mm_store_ps (tab,x);
    for (i=0;i<4;i++)
        printf ("%8.2f ",tab[i]);
}

#endif


#ifdef ALGEBRA_AVX

void _mm256_print_ps(__m256 x) {
    size_t i;
    float tab[8];
    _mm256_store_ps (tab,x);
    for (i=0;i<8;i++)
        printf ("%8.2f ",tab[i]);
}

#endif

void	algebra_show_arch()
{
	#ifdef ALGEBRA_SSE2
	printf("sse2\t: yes\n");
	#else
	printf("sse2\t: no\n");
	#endif
	#ifdef ALGEBRA_SSSE3
	printf("ssse3\t: yes\n");
	#else
	printf("ssse3\t: no\n");
	#endif
	#ifdef ALGEBRA_SSE4
	printf("sse4\t: yes\n");
	#else
	printf("sse4\t: no\n");
	#endif
	#ifdef ALGEBRA_AVX
	printf("avx\t: yes\n");
	#else
	printf("avx\t: no\n");
	#endif
}

