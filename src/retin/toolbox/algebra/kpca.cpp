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
 * \file KPCA.cpp
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#include "kpca.h"

#include "vector_float.h"

#include <algorithm>

using namespace std;

namespace algebra {

KPCA::KPCA (size_t dim,double* gram,bool verbose) :
    dim(dim),gram(gram),verbose(verbose),featureDim(0)
{

}

KPCA::~KPCA() {
    
}

void    KPCA::computeEigens(bool centering) {
    eigens.clear();
    std::vector<double> vectors(dim*dim),lambdas(dim);

    memcpy(&vectors[0],gram,sizeof(double)*dim*dim);
    if (centering)
    {
        std::cout << "kpca => centering" << std::endl;
        matrix_sym_centering_double (&vectors[0],dim,NULL);
    }

    matrix_eigSym_double(&vectors[0],&lambdas[0],dim);
    
    double max = 0;
    for (size_t r=0;r<dim;r++) {
        if (fabs(lambdas[r]) > max)
            max = fabs(lambdas[r]);
    }
    if (max < 1E-7)
        return;
    for (size_t r=0;r<dim;r++) {
        if (fabs(lambdas[r]) > 1E-6*max) {
            eigens.push_back(Eigen(lambdas[r],&vectors[r*dim],dim));
        }
    }
}

void    KPCA::sort(size_t maxDim,double maxEnergy) {
    if (maxDim == 0)
        maxDim = dim;

    double sum = 0;
    for (size_t r=0;r<eigens.size();r++) {
        sum += fabs(eigens[r].getValue());
    }
    fullEnergy = sum;
    std::sort(eigens.begin(),eigens.end());

    size_t dim = 0;
    double sumper = 0;
    for (size_t r=0;r<eigens.size();r++) {
        double per = 100*fabs(eigens[r].getValue())/sum;
        if (per < 1E-7)
            break;
        sumper += per;
        if (verbose)
            std::cout << per << "% ";
        dim ++;
        if (dim >= maxDim || sumper >= maxEnergy)
            break;
    }
    if (verbose) {
        std::cout << "= " << sumper << "% <= " << maxEnergy << "%" << std::endl;
        std::cout << "dim = " << dim << std::endl;
    }
    eigens.resize(dim);
}

void    KPCA::run(size_t maxDim,double maxEnergy,bool centering) {
    computeEigens(centering);
    sort(maxDim,maxEnergy);
}

void    KPCA::setEigens (const double* values,const double* vectors,size_t count) {
    eigens.clear();
    for (size_t i=0;i<count;i++) {
        eigens.push_back(Eigen(values[i],vectors+i*dim,dim));
    }
}


void    KPCA::compPrimalBegin(size_t featureDim_,size_t firstEigen_,size_t eigenCount_)
{
    featureCount = 0;
    featureDim = featureDim_;
    firstEigen = firstEigen_;
    eigenCount = eigenCount_;
    if (firstEigen+eigenCount > eigens.size())
        retinThrowException("Invalid eigen dims");
    primalFeature.resize(featureDim);
    primalMean.resize(featureDim);
    primalValues.resize(eigenCount);
    primalVectors.resize(featureDim*eigenCount);

}

void    KPCA::compPrimalAdd(const float* feature,size_t featureDim_,size_t firstEigen_,size_t eigenCount_)
{
    if (featureDim == 0)
        compPrimalBegin(featureDim_,firstEigen_,eigenCount_);
    if (featureDim != featureDim_)
        retinThrowException2("Bad featureDim %d != %d",featureDim_,featureDim);
    if (firstEigen != firstEigen_)
        retinThrowException2("Bad firstEigen %d != %d",firstEigen_,firstEigen);
    if (eigenCount != eigenCount_)
        retinThrowException2("Bad eigenCount %d != %d",eigenCount_,eigenCount);
    if (featureCount > dim)
        retinThrowException1("Wrong number of features (%ld expected)",dim);
    vector_add_float(&primalMean[0],feature,featureDim_);
    for (size_t j=0;j<eigenCount;j++) {
        double eig = 1.0 / sqrt(eigens[firstEigen+j].getValue());
        const double* eigv = eigens[firstEigen+j].getVector();
        vector_addm_float(&primalVectors[j*featureDim],eigv[featureCount] * eig,feature,featureDim);
    }
    featureCount ++;
}

void    KPCA::compPrimalEnd ()
{
    if (featureCount != dim)
        retinThrowException("Wrong number of features");
    for (size_t i=0;i<featureDim;i++)
        primalMean[i] /= dim;
    for (size_t j=0;j<eigenCount;j++)
        primalValues[j] = (float)(eigens[firstEigen+j].getValue() / dim);
}


void    KPCA::compPrimal(const float* features,size_t featureDim_,size_t featureCount_)
{
    if (featureCount_ != dim)
        retinThrowException2("Bad featureCount %d != %d",featureCount_,dim);

    compPrimalBegin(featureDim_,0,eigens.size());
    for (size_t i=0;i<featureCount_;i++) {
        compPrimalAdd(features+i*featureDim_,featureDim_,0,eigens.size());
    }
    compPrimalEnd();

    /*size_t p = eigens.size();
    vector<double> ls(p);
    vector<const double*> vs(p);
    for (size_t j=0;j<eigens.size();j++) {
        ls[j] = 1.0 / sqrt(eigens[j].getValue());
        vs[j] = eigens[j].getVector();
        values[j] = eigens[j].getValue() / dim;
    }
    for (size_t k=0;k<dim;k++) {
        for (size_t i=0;i<featureDim;i++) {
            double x = features[i+k*featureDim];
            for (size_t j=0;j<p;j++) {
                temp[i+j*featureDim] += x * vs[j][k] * ls[j];
            }
        }
    }
    for (size_t i=0;i<temp.size();i++)
        vectors[i] = (float)temp[i];*/

    /*for (size_t j=0;j<eigens.size();j++) {
        double l = 1.0 / sqrt(eigens[j].getValue());
        const double* v = eigens[j].getVector();
        values[j] = eigens[j].getValue() / dim;
        for (size_t i=0;i<featureDim;i++) {
            double w = 0;
            for (size_t k=0;k<dim;k++) {
                w += features[i+k*featureDim] * v[k] * l;
            }
            vectors[i+j*featureDim] = w;
        }
    }*/
}



}
