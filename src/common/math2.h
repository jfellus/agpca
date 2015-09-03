/*
 * math.h
 *
 *  Created on: 8 nov. 2013
 *      Author: jfellus
 */

#ifndef GOSSIP_MATH_H_
#define GOSSIP_MATH_H_


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "matrix.h"

float frand(float min, float max) {
	return ((float)rand())/INT_MAX*(max-min)+min;
}

void randvec(float* v, int n, float min, float max) {
	for(int i=0; i<n; i++) {
		v[i] = frand(min,max);
	}
}

void randvec_gaussian(float* v, int d, float* mu, float* sigma) {

}

#endif /* GOSSIP_MATH_H_ */
