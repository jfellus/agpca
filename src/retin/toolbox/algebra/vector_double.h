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

#ifndef __algebra_vector_double_h__
#define __algebra_vector_double_h__

#include "core.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Calcule la norme L1. */
double	vector_n1_double (const double* v,size_t n);
/*! Calcule la norme L2. */
double	vector_n2_double (const double* v,size_t n);
/*! Calcule la norme L2 au carré. */
double	vector_n2p2_double (const double* v,size_t n);
/*! Calcule la distance L2 entre deux vectors. */
double	vector_l2_double (const double* v1,const double* v2,size_t n);
/*! Calcule la distance L2 au carré entre deux vectors. */
double	vector_l2p2_double (const double* v1,const double* v2,size_t n);

/*! Calcule le produit scalaire. */
double	vector_ps_double (const double* v1,const double* v2,size_t n);
/*! Calcule le produit scalaire. version bench */
double	vector_ps_double_basic (const double* v1,const double* v2,size_t n);

/*! Calcule la valeur max. */
double	vector_max_double (const double* v,size_t n);

/*! Calcule v1 += v2. */
void	vector_add_double (double* v1,const double* v2,size_t n);
/*! Calcule v1 += lambda*v2. */
void	vector_addm_double (double* v1,double lambda,const double* v2,size_t n);
void	vector_addm_double_basic (double* v1,double lambda,const double* v2,size_t n);
/*! Calcule v1 -= v2. */
void	vector_sub_double (double* v1,const double* v2,size_t n);
/*! Calcule v1 *= v2. */
void	vector_mul_double (double* v1,double* v2,size_t n);
/*! Calcule v1 /= v2. */
void	vector_div_double (double* v1,double* v2,size_t n);

/*! Calcule v = v1+alpha*v2. */
void	vector_linear_double (double* v,const double* v1,const double alpha,const double* v2,size_t n);
/*! Calcule mean += lambda*v2 et var += lambda*v2^2 (pour calculer moyenne et variance par la suite). */
void	vector_add_stats_double (double* mean,double* var,double lambda,const float* v,size_t n);

/*! Calcule v = r. */
void	vector_scpy_double (double* v,double r,size_t n);
/*! Calcule v += r. */
void	vector_sadd_double (double* v,double r,size_t n);
/*! Calcule v -= r. */
void	vector_ssub_double (double* v,double r,size_t n);
/*! Calcule v *= r. */
void	vector_smul_double (double* v,double r,size_t n);
/*! Calcule v /= r. */
void	vector_sdiv_double (double* v,double r,size_t n);


/*! Calcule p[i] = exp(p[i]) / sum_j(exp(p[j]) de manière prudente, pour éviter les erreurs de précision */
void    vector_exp_proba_double (double* p,size_t n);


/*! Calcule la valeur min. */
double	vector_min_double (const double* v,size_t n);
/*! Cacule la somme des valeurs. */
double	vector_sum_double (const double* v,size_t n);

#ifdef __cplusplus
}
#endif

#endif
