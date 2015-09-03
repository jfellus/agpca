/*
 * math.h
 *
 *  Created on: 8 nov. 2013
 *      Author: jfellus
 */

#ifndef GOSSIPS_MATH_H_
#define GOSSIPS_MATH_H_


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

class MatrixFloat;
class MatrixDouble;


float frand();
float frand(float min, float max);

float rand_gaussian(float mu, float sigma);
float rand_exp(float lambda);

void randvec(float* v, int n, float min, float max);
void randvec(double* v, int n, double min, double max);

void randvec_gaussian(float* v, int d);
void randvec_gaussian(float* v, MatrixFloat& mu, MatrixFloat& sigma);
void rand_covariance(MatrixFloat& cov, float s, int rank = 5);

void randvec_gaussian(double* v, int d);
void randvec_gaussian(double* v, MatrixDouble& mu, MatrixDouble& sigma);
void rand_covariance(MatrixDouble& cov, double s, int rank = 5);

#endif /* GOSSIPS_MATH_H_ */
