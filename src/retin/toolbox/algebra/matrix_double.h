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

#ifndef __algebra_matrix_double_h__
#define __algebra_matrix_double_h__

#include "core.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! Calcule C += A*B. */
void	matrix_CpAB_double (double* C,double* A,double* B,size_t n,size_t p,size_t m);
/*! Calcule C += A'*B. */
void	matrix_CpAtB_double (double* C,double* A,double* B,size_t n,size_t p,size_t m);
/*! Calcule C += A*B'. */
void	matrix_CpABt_double (double* C,double* A,double* B,size_t n,size_t p,size_t m);
/*! Calcule C += a*a' (pratique pour matrice covariance). */
void	matrix_Cpaat_double (double* C,const double* a,size_t n);
/*! Calcule C += A*A'. */
/*! C de taille (n,n)
 *  A de taille (n,p)
 */
void	matrix_CpAAt_double (double* C,const double* A,size_t n,size_t p);
/*! Calcule C += A*B*A'. */
/*! r est un vecteur temporaire de taille m. */
void	matrix_CpABAt_double (double* C,double* A,double* B,size_t n,size_t m,double* r);
/*! Calcule C += A'*B*A. */
/*! r est un vecteur temporaire de taille m. */
void	matrix_CpAtBA_double (double* C,double* A,double* B,size_t n,size_t m,double* r);
/*! Calcule C += A*d*A' (d vecteur diagonale). */
void	matrix_CpAdAt_double (double* C,double* A,double* d,size_t n,size_t m);


/*! Calcule Y = X'. */
void	matrix_t_double (double* Y,double* X,size_t n,size_t m);
/*! Cacule s_{1j} = sum_i X_{ij}. */
void	matrix_sum_double (double* s,double* X,size_t n,size_t m);
/*! Cacule s_{i} = sum_j X_{ij}. */
void	matrix_sumt_double (double* s,double* X,size_t n,size_t m);
/*! Cacule s_{i} = 1/m sum_i X_{ij}. */
void	matrix_meant_double(double* s,const double* X,size_t n,size_t m);
/*! Calcule la distance L2 entre deux matrices. */
double	matrix_l2_double (const double* A,const double* B,size_t n,size_t m);
/*! Calcule la distance L2 au carré entre deux matrices. */
double	matrix_l2p2_double (const double* A,const double* B,size_t n,size_t m);

/*! Décomposition LU. */
int	matrix_lu_double (double* M,size_t* idx,double* d,size_t n);
/*! Backsubstition LU. */
void	matrix_lub_double (double* b,double* M,size_t *idx,size_t n);
/*! Inversion. */
/*! Y: matrix(n,n) inverse.\n
    M: matrix(n,n) à inverser (modifiée).*/
int	matrix_inv_double (double* Y,double* M,size_t n);
/*! Inversion et multiplication Y = inv(M)*B. */
/*! Y: matrix(n,n) inverse.\n
    M: matrix(n,n) à inverser (modifiée).\n
    B: matrix(n,m) à multiplier.*/
int	matrix_invm_double (double* Y,double* M,double* B,size_t n,size_t m);
/*! Décomposition Cholesky. */
size_t	matrix_ll_double (double* L,double* A,size_t n,size_t m);
/*! Décomposition QR. */
size_t	matrix_qr_double (double* A,double* R,size_t n,size_t m);
/*! Décomposition QT. */
void	matrix_qtq_double (double* A,double* d,double* e,size_t n);
/*! Décomposition QT partielle (que le T). */
void	matrix_qtq_t_double (double* A,double* d,double* e,size_t n);
/*! Décomposition complète en valeur et vecteurs propres d'une matrix trigonale. */
void	matrix_eigTrig_double (double* A,double* d,double* e,size_t n);
/*! Calcul des valeurs propres d'une matrix trigonale. */
void	matrix_eigTrig_values_double (double* d,double* e,size_t n);
/*! Décomposition complète en valeur et vecteurs propres d'une matrix symétrique. */
void	matrix_eigSym_double (double* A,double* d,size_t n);
/*! Calcul de la plus grande valeur propre. */
int	matrix_eigMax_double (double* l,double* v,double* M,size_t n);
/*! Centre une matrice symétrique. */
/*! buf: tampon de n doubles (alloué en interne si NULL)
 */
void	matrix_sym_centering_double (double* M,size_t n,const double* buf);

/*! Sauvegarde le contenu d'une matrice dans un fichier */
int     matrix_save_double(const char* file_name,const double* M,size_t n,size_t m);

/*! Recherche du vecteur le plus proche dans un dictionnaire. */
/*!
 * sample : vecteur à classifier
 * dict : matrice(n,m) dictionnaire
 * n: taille des vecteurs
 * m: nombre de vecteurs dans le dico.
 */
size_t matrix_argmin_l2_double (const double* sample,const double* dict,size_t n,size_t m);


/*! Trie dans l'ordre décroissant les valeurs propres. */
void	matrix_sortEig_double (double* A,double* d,size_t n);


/*! Calcule la norme L2. */
double	matrix_n2_double (const double* v,size_t n,size_t m);

#ifdef __cplusplus
}
#endif

#endif
