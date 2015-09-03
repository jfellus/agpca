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

#ifndef __algebra_vector_uchar_h__
#define __algebra_vector_uchar_h__

#include "core.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Calcule le produit scalaire. */
unsigned int	vector_ps_uchar (const unsigned char* v1,const unsigned char* v2,size_t n);
/*! Calcule la distance L1 entre deux vectors. */
unsigned int	vector_l1_uchar (const unsigned char* v1,const unsigned char* v2,size_t n);
/*! Calcule la distance L2 au carré entre deux vectors. */
unsigned int	vector_l2p2_uchar (const unsigned char* v1,const unsigned char* v2,size_t n);


#ifdef __cplusplus
}
#endif

#endif
