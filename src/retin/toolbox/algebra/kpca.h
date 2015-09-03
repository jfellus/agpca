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
 * \file KPCA.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __algebra_KPCA_h__
#define __algebra_PCA_h__

#include "vector_double.h"
#include "matrix_double.h"

#include "retin/toolbox/core/auto_array_ptr.h"


namespace algebra {

class KPCA {

public:
    class Eigen
    {
        double                          l;
        std::vector<double>	v;
    public:
        Eigen(double l=0,const double* pv=NULL,size_t n=0) : l(l),v(n) {
                for (size_t i=0;i<n;i++) {
                    v[i] = pv[i];
                }
        }
        bool	operator<(const Eigen& e) const {
                return fabs(l) > fabs(e.l);
        }
        double          getValue() const { return l; }
        const double*   getVector() const { return &v[0]; }
    };

    // Ne prends pas possession du buffer "gram" (il faut le libérer soit même).
    KPCA (size_t dim,double* gram,bool verbose=false);
    virtual ~KPCA();
    void                setVerbose(bool b) { verbose = b; }
   

    void                run (size_t maxDim,double maxEnergy,bool centering);

    size_t                          getDim() const { return dim; }
    const double*                   getGram() const { return gram; }
    const std::vector<Eigen>&       getEigens() const { return eigens; }
    double                          getEnergy() const { return fullEnergy; }

    void                setEigens (const double* values,const double* vectors,size_t count);

    //! Calcule toutes les valeurs propres. A utiliser à la place de "run".
    void                computeEigens(bool centering);
    void                sort(size_t maxDim,double maxEnergy=100);



protected:
    size_t  dim;
    double fullEnergy;
    double* gram;
    bool verbose;
    std::vector<Eigen> eigens;

// Calculs du primal en streaming
// ATTENTION IL FAUT ABSOLUMENT DONNER LES FEATURES DANS LE MEME ORDRE
// QUE LA MATRICE DE GRAM !!!
protected:
    size_t  featureDim,featureCount;
    size_t  firstEigen,eigenCount;
    retin::auto_array_ptr<float> primalFeature;
    retin::auto_array_ptr<float> primalMean;
    retin::auto_array_ptr<float> primalValues;
    retin::auto_array_ptr<float> primalVectors;
public:
    // Calcul en streaming (feature par feature)
    // Ne fera les calculs que pour les valeurs propres entre firstEigen et firstEigen+eigenCount-1
    void        compPrimalBegin (size_t featureDim_,size_t firstEigen_,size_t eigenCount_);
    void        compPrimalAdd (const float* feature,size_t featureDim_,size_t firstEigen_,size_t eigenCount_);
    void        compPrimalEnd ();

    size_t      getPrimalDim() const { return featureDim; }
    size_t      getFirstSelectedEigen() const { return firstEigen; }
    size_t      getSelectedEigenCount() const { return eigenCount; }
    // Donne le buffer avec la moyenne (attention il faudra le libérer soit même)
    float*      releasePrimalMean() { return primalMean.release(); }
    // Donne le buffer avec les valeurs propres (attention il faudra le libérer soit même)
    float*      releasePrimalValues() { return primalValues.release(); }
    // Donne le buffer avec les vecteurs propres (attention il faudra le libérer soit même)
    float*      releasePrimalVectors() { return primalVectors.release(); }

    // Calcul en une passe (appelle les autres)
    void        compPrimal(const float* features,size_t featureDim,size_t featureCount);
};


}


#endif
