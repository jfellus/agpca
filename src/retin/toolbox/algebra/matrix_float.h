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

#ifndef __algebra_matrix_float_h__
#define __algebra_matrix_float_h__

#include "vector_float.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Calcule C += A*B. */
/*! C de taille (n,m)
 *  A de taille (n,p)
 *  B de taille (p,m)
 */
void	matrix_CpAB_float (float* C,const float* A,const float* B,size_t n,size_t p,size_t m);
/*! Calcule C += A'*B. */
void	matrix_CpAtB_float (float* C,const float* A,const float* B,size_t n,size_t p,size_t m);
/*! Calcule C += A*B'. */
/*! C de taille (n,m)
 *  A de taille (p,n)
 *  B de taille (p,m)
 */
void	matrix_CpABt_float (float* C,const float* A,const float* B,size_t n,size_t p,size_t m);
/*! Calcule C += A*A'. */
/*! C de taille (n,n)
 *  A de taille (n,p)
 */
void	matrix_CpAAt_float (float* C,const float* A,size_t n,size_t p);

/*! Calcule C += a*a' (pratique pour matrice covariance). */
void	matrix_Cpaat_float (float* C,const float* a,size_t n);

/*! Calcule C += A*B*A'. */
/*! r est un vecteur temporaire de taille m. */
void	matrix_CpABAt_float (float* C,const float* A,const float* B,size_t n,size_t m,float* r);
/*! Calcule C += A'*B*A. */
/*! r est un vecteur temporaire de taille m. */
void	matrix_CpAtBA_float (float* C,const float* A,const float* B,size_t n,size_t m,float* r);
/*! Calcule C += A*d*A' (d vecteur diagonale). */
void	matrix_CpAdAt_float (float* C,const float* A,const float* d,size_t n,size_t m);

/*! Calcule Y = X'. */
void	matrix_t_float (float* Y,const float* X,size_t n,size_t m);
/*! Calcule X = X' */
void    matrix_selft_float(float* X,size_t n,size_t m);
/*! Calcule Y = X(i:i+n,j:j+m). */
void	matrix_cpy_float (float* Y,size_t iy,size_t jy,const float* X,size_t ix,size_t jx,size_t n,size_t m,size_t nx,size_t ny);
/*! Calcule Y(i,j) = X((i+di)%n,(j+dj)%m). */
void	matrix_rot_float (float* Y,const float* X,size_t di,size_t dj,size_t n,size_t m);
/*! Cacule s_{1j} = sum_i X_{ij}. */
void	matrix_sum_float (float* s,const float* X,size_t n,size_t m);
/*! Cacule s_{i} = sum_j X_{ij}. */
void	matrix_sumt_float (float* s,const float* X,size_t n,size_t m);
/*! Cacule s_{i} = 1/m sum_i X_{ij}. */
void	matrix_meant_float (float* s,const float* X,size_t n,size_t m);
/*! Calcule le produit scalaire (frobenuis) entre deux matrices. */
double	matrix_ps_float (const float* A,const float* B,size_t n,size_t m);
/*! Calcule le produit scalaire (frobenuis) entre xx' et yy'. */
double	matrix_ps_xxt_yyt_float (const float* x,const float* y,size_t n);
/*! Calcule le produit scalaire (frobenuis) entre XX' et YY'. */
double	matrix_ps_XXt_YYt_float (const float* X,size_t p,const float* Y,size_t q,size_t n);
/*! Calcule la norme L2. */
double	matrix_n2_float (const float* v,size_t n,size_t m);
/*! Calcule la distance L2 entre deux matrices. */
double	matrix_l2_float (const float* A,const float* B,size_t n,size_t m);
/*! Calcule la distance L2 au carré entre deux matrices. */
double	matrix_l2p2_float (const float* A,const float* B,size_t n,size_t m);

/* Calculs centrés, avec H = Identité - 1/n */

/*! Calcule la norme centré tr(HKHK). */
double  matrix_cn2_float (const float* K,size_t n);
/*! Calcule le produit scalaire centré tr(HKHL). */
/*! K et L sont deux matrices de taille n x n */
double  matrix_cps_float(const float* K,const float* L,size_t n);
/*! Calcule le produit scalaire centré entre xx' et yy'. */
/*! x et y sont des vecteurs de taille n */
double	matrix_cps_xxt_yyt_float (const float* x,const float* y,size_t n);
/*! Calcule le produit scalaire centré entre XX' et yy'. */
/*! X est une matrice de p vecteurs de taille n
 *  y est un vecteur de taille n */
double	matrix_cps_XXt_yyt_float (const float* X,size_t p,const float* y,size_t n);
/*! Calcule le produit scalaire centré entre XX' et YY'. */
/*! X est une matrice de p vecteurs de taille n
 *  Y est une matrice de q vecteurs de taille n */
double	matrix_cps_XXt_YYt_float (const float* X,size_t p,const float* Y,size_t q,size_t n);
/**
 * Calcul la matrice diagonal: diag(HXX'Hyy')
 * @param X est une matrice de p vecteurs de taille n
 * @param y est un vecteur de taille n
 * @param res est le vecteur resultat de taille n
 */
void    matrix_diag_XXt_yyt_float (const float* X,size_t p,const float* y,size_t n,float *res);

/*! Calcule le produit scalaire centré entre X'X et Y'Y. */
/*! X est une matrice de n vecteurs de taille p
 *  Y est une matrice de n vecteurs de taille q */
double	matrix_cps_XtX_YtY_float (const float* X,size_t p,const float* Y,size_t q,size_t n);

/*! Décomposition LU. */
int     matrix_lu_float (float* M,size_t* idx,float* d,size_t n);
/*! Backsubstition LU. */
void	matrix_lub_float (float* b,float* M,size_t *idx,size_t n);
/*! Inversion. */
/*! Y: matrix(n,n) inverse.\n
    M: matrix(n,n) à inverser (modifiée).*/
int     matrix_inv_float (float* Y,float* M,size_t n);
/*! Inversion et multiplication Y = inv(M)*B. */
/*! Y: matrix(n,n) inverse.\n
    M: matrix(n,n) à inverser (modifiée).\n
    B: matrix(n,m) à multiplier.*/
int     matrix_invm_float (float* Y,float* M,float* B,size_t n,size_t m);
/*! Inverse de la matrice de covariance C_{ij} = X_{i}^t X^{j} (=> Moindre carrés/LLE). */
/*! invC : matrix(m,m) inverse
 *  X : matrix(n,m) de m vecteurs de dimension n (modifiée apres appel)
 *  Renvoie zero si n<m, sinon le rang de la matrice de covariance */
size_t	matrix_inv_cov_float (float* invC,float* X,size_t n,size_t m);
/*! Approximation LLE. */
/*! y : vecteur(m) coefficients
 *  x : vecteur(m) à approximer
 *  X : matrix(n,m) base de m vecteurs de dimension n
 *  invC : matrix(m,m) inverse des covariances de X */
void	matrix_lle_float(float* y,const float* x,const float* X,const float* invC,size_t n,size_t m);
/*! Décomposition Cholesky. */
size_t	matrix_ll_float (float* L,float* A,size_t n,size_t m);
/*! Inversion d'une matrice triangulaire inférieure. */
/*! Y: matrix(n,n) inverse.\n
    L: matrix(n,n) à inverser.*/
void	matrix_inv_r_float (float* Y,float* L,size_t n);
/*! Décomposition QR. */
size_t	matrix_qr_float (float* A,float* R,size_t n,size_t m);
/*! Décomposition QT. */
void	matrix_qtq_float (float* A,float* d,float* e,size_t n);
/*! Décomposition QT partielle (que le T). */
void	matrix_qtq_t_float (float* A,float* d,float* e,size_t n);
/*! Décomposition complète en valeur et vecteurs propres d'une matrix trigonale. */
void	matrix_eigTrig_float (float* A,float* d,float* e,size_t n);
/*! Calcul des valeurs propres d'une matrix trigonale. */
void	matrix_eigTrig_values_float (float* d,float* e,size_t n);
/*! Décomposition complète en valeur et vecteurs propres d'une matrix symétrique. */
void	matrix_eigSym_float (float* A,float* d,size_t n);
/*! Trie dans l'ordre décroissant les valeurs propres. */
void	matrix_sortEig_float (float* A,float* d,size_t n);
/*! Calcul de la plus grande valeur propre. */
int     matrix_eigMax_float (float* l,float* v,float* M,size_t n);
int     matrix_eigMaxSym_float (float* l,float* v,float* M,size_t n);
/*! Centre une matrice symétrique. */
/*! buf: tampon de n floats (alloué en interne si NULL)
 */
void	matrix_sym_centering_float (float* M,size_t n,const float* buf);


/*! Recherche du vecteur le plus proche dans un dictionnaire. */
/*!
 * sample : vecteur à classifier
 * dict : matrice(n,m) dictionnaire
 * n: taille des vecteurs
 * m: nombre de vecteurs dans le dico.
 */
size_t  matrix_argmin_l2_float (const float* sample,const float* dict,size_t n,size_t m);

#ifdef __cplusplus
}
#endif

#endif
