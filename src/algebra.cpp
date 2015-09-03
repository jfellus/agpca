/*
 * algebra.cpp
 *
 *  Created on: 13 nov. 2013
 *      Author: jfellus
 */

#include <math.h>



int discard_after(float* U,float* L, int NB_PROJS) {
	//DBG("DISCARDDDD!!!! " << NB_PROJS);
	for(int i=NB_PROJS; i<D; i++) {
		for(int j=0; j<D; j++) {
			U[i*D + j] = L[i] = 0;
		}
	}
	return NB_PROJS;
}

int discard_after_energy(float* U,float* L, float keep_energy) {
	DBG("DISCARDDDD ENERGY!!!! " << q);
	float sum = vector_sum_float(L,D);
	float e = 0;
	int nb_projs_kept = D;
	for(int i=0; i<D; i++) {
		if(e >= keep_energy) {
			if(nb_projs_kept==D) nb_projs_kept = i;
			for(int j=0; j<D; j++) U[i*D+j] = 0;
			L[i] = 0;
		}
		e += (L[i])/sum;
	}
	return nb_projs_kept;
}


