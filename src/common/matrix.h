/*
 * matrix.h
 *
 *  Created on: 23 janv. 2014
 *      Author: jfellus
 */

#ifndef MATRIX_H_
#define MATRIX_H_


#ifdef USE_MATRIX_DOUBLE
	#include "matrix_double.h"
	typedef MatrixDouble Matrix;
#else
	#include "matrix_float.h"
	typedef MatrixFloat Matrix;
#endif

#endif /* MATRIX_H_ */
