/*
 * the_eigen.h
 *
 *  Created on: 4 sept. 2015
 *      Author: jfellus
 */

#ifndef THE_EIGEN_H_
#define THE_EIGEN_H_


void QR(float* A, float* R, size_t n, size_t m);
void QR(float* A, float* Q, float* R, size_t n, size_t m) ;
void qr(float* A, float* Q, float* r, size_t n, size_t m);
void qr(float* A, float * r, size_t n, size_t m);


void eigendecompose(float* A, float* U, float* L, size_t n);



#endif /* THE_EIGEN_H_ */
