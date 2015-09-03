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
 * \file matrix.h
 * \author Philippe H. Gosselin
 * \version 1.0
 */

#ifndef __algebra_matrix_hpp__
#define __algebra_matrix_hpp__

#include "matrix_float.h"
#include "matrix_double.h"

template<typename U,typename T> void matrix_convert(U* v1,const T* v2,size_t n,size_t m) { vector_convert(v1,v2,n*m); }


#include <iostream>

template<typename T>
void matrix_print(std::ostream& out,const T* M,size_t n,size_t m,const char* format="%8.2f ")
{
    size_t i,j;
    char tmp[0x100];
    for (i=0;i<n;i++) {
        for (j=0;j<m;j++) {
            sprintf(tmp,format,M[i+j*n]);
            out << tmp;
        }
        out << std::endl;
    }
}


#endif
