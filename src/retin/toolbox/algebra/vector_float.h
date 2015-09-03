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

#ifndef __algebra_vector_float_h__
#define __algebra_vector_float_h__

#include "core.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Matrice aléatoire uniforme */
void    vector_rand_float (float* v,float fMin,float fMax,size_t n);

/*! Calcule le produit scalaire. */
float	vector_ps_float (const float* v1,const float* v2,size_t n);
/*! Calcule le produit scalaire. */
/* Version non optimizée, pour les benchmarks */
float	vector_ps_float_basic (const float* v1,const float* v2,size_t n);
/*! Calcule le produit scalaire par parquet de 32 valeurs.
 * Les deux vecteurs v1 et v2 sont de taille 32*q
 * q doit être non nul
 * Attention : les données doivent être alignées */
float	vector_ps_float_aligned_32 (const float* v1,const float* v2,size_t q);

/*! Cacule la somme des valeurs. */
float	vector_sum_float (const float* v,size_t n);
/*! Cacule la somme des valeurs. */
/* Version non optimizée, pour les benchmarks */
float	vector_sum_float_basic (const float* v,size_t n);

/*! Calcule la norme L1. */
float	vector_n1_float (const float* v,size_t n);
/*! norm L1, version basic pour benchmark */
float	vector_n1_float_basic (const float* v,size_t n);
/* norme L1, version vectorisée */
float vector_n1_float_aligned_32(const float * v, size_t q);

/*! Calcule la norme L2. */
float	vector_n2_float (const float* v,size_t n);
/*! Calcule la norme L2 au carré. */
float	vector_n2p2_float (const float* v,size_t n);
/*! Calcule la norme L2 au carré.version basic pour bench */
float	vector_n2p2_float_basic (const float* v,size_t n);

/*! Calcule la distance L1 entre deux vectors. */
float	vector_l1_float (float* v1,float* v2,size_t n);
/*! Calcule la distance L2 entre deux vectors. */

float	vector_l2_float (const float* v1,const float* v2,size_t n);
/*! Calcule la distance L2 au carré entre deux vectors. */
float	vector_l2p2_float (const float* v1,const float* v2,size_t n);
/*! Calcule la distance L2 au carré entre deux vectors. version basic pour bench */
float	vector_l2p2_float_basic (const float* v1,const float* v2,size_t n);

/*! Calcule la distance Chi1 entre deux vectors. */
float	vector_chi1_float (float* v1,float* v2,size_t n);
/*! Calcule la distance Chi1 entre deux vectors translatés de t. */
float	vector_chi1t_float (float* v1,float* v2,float* t,size_t n);
/*! Calcule la distance Chi2 entre deux vectors. */
float	vector_chi2_float (const float* v1,const float* v2,size_t n);
/*! Calcule la similarité cos entre deux vectors. */
float	vector_cos_float (const float* v1,const float* v2,size_t n);
/*! Calcule la distance Linf entre deux vectors. */
float	vector_linf_float (float* v1,float* v2,size_t n);
/*! Calcule la valeur max. */
float	vector_max_float (float* v,size_t n);
/*! Calcule la valeur min. */
float	vector_min_float (float* v,size_t n);
/*! Calcule l'écart type. */
float	vector_std_float (const float* v1,const float* v2,size_t n);
/*! Calcule l'indice de la valeur max. */
size_t	vector_argmax_float (float* v,size_t n);
/*! Calcule l'indice de la valeur min. */
size_t	vector_argmin_float (float* v,size_t n);

/*! Calcule v[i] = abs(v[i]). */
void	vector_abs_float (float* v,size_t n);
/*! Calcule v[i] = sgn(v[i])*sqrt(abs(v[i]),p). */
void	vector_sqrt_float (float* v,size_t n);
/*! Calcule v[i] = sgn(v[i])*powf(abs(v[i]),p). */
void	vector_pow_float (float* v,float p,size_t n);
/*! Règle de 3 pour avoir des valeurs entre rMin et rMax. */
void	vector_rescale_float (float* v,float rMin,float rMax,size_t n);
/*! Calcule la norme complexe z[i] = sqrt(x[i]*x[i]+y[i]*y[i]). */
void	vector_cn2_float (float* z,float* x,float* y,size_t n);
/*! Calcule la norme quaternionique. */
void    vector_qn2_float (float* z,float* x,float* yi,float* yj,float* yk,size_t n);

/*! Calcule v1 = lambda*v2. */
void	vector_cpym_float (float* v1,float lambda,const float* v2,size_t n);
/*! Calcule v1 += v2. */
void	vector_add_float (float* v1,const float* v2,size_t n);
/*! Calcule v1 += lambda*v2. */
void	vector_addm_float (float* v1,float lambda,const float* v2,size_t n);
void	vector_addm_float_basic (float* v1,float lambda,const float* v2,size_t n);
/*! Calcule v1 -= v2. */
void	vector_sub_float (float* v1,const float* v2,size_t n);
/*! Calcule v1 -= v2. Version non optimisée pour les benchmarks */
void	vector_sub_float_basic (float* v1,const float* v2,size_t n);
/*! Calcule v1 *= v2. */
void	vector_mul_float (float* v1,float* v2,size_t n);
/*! Calcule v1 /= v2. */
void	vector_div_float (float* v1,float* v2,size_t n);

/*! Calcule v = v1+alpha*v2. */
void	vector_linear_float (float* v,const float* v1,const float alpha,const float* v2,size_t n);

/*! Calcule v = r. */
void	vector_scpy_float (float* v,float r,size_t n);
/*! Calcule v += r. */
void	vector_sadd_float (float* v,float r,size_t n);
/*! Calcule v -= r. */
void	vector_ssub_float (float* v,float r,size_t n);
/*! Calcule v *= r. */
void	vector_smul_float (float* v,float r,size_t n);
/*! Calcule v /= r. */
void	vector_sdiv_float (float* v,float r,size_t n);


/*! Vérifie que les valeurs du vecteur sont entières. */
int     vector_isfinite_float (const float* v,size_t n);

size_t  vector_argmin_l2_float (const float* x,const float* C,size_t n,size_t q,const float* normC);


#ifdef __cplusplus
}
#endif

#endif
