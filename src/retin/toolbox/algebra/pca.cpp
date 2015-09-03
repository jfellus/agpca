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
 * \file PCA.cpp
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#include "pca.h"

#include<algorithm>

using namespace std;

namespace algebra {

void    PCA::center() {
    if (count == 0)
        return;
    vector_sdiv_double(&mean[0],count,dim);
    for (size_t s=0;s<dim;s++) {
        for (size_t r=0;r<dim;r++) {
            covs[r+s*dim] = covs[r+s*dim]/count - mean[r]*mean[s];
        }
    }
}

void    PCA::computeEigens() {
    eigens.clear();
    if (count == 0)
        return;
    std::vector<double> vectors(dim*dim),lambdas(dim);
    memcpy(&vectors[0],&covs[0],dim*dim*sizeof(double));
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

void    PCA::run(size_t maxDim,double maxEnergy) {
    if (verbose)
        std::cout << "Calcul des projections..." << std::endl;
    computeEigens();
    sort(maxDim,maxEnergy);
}

void    PCA::sort(size_t maxDim,double maxEnergy) {
    if (maxDim == 0)
        maxDim = dim;

    double sum = 0;
    for (size_t r=0;r<eigens.size();r++) {
        sum += fabs(eigens[r].getValue());
    }
    totEnergy = sum;
    std::sort(eigens.begin(),eigens.end());

    size_t dim = 0;
    double sumper = 0;
    for (size_t r=0;r<eigens.size();r++) {
        double per = 100*fabs(eigens[r].getValue())/sum;
        if (per < 1E-7)
            break;
        sumper += per;
        if (verbose)
            std::cout << per << "%% ";
        dim ++;
        if (dim >= maxDim || sumper >= maxEnergy)
            break;
    }
    if (verbose) {
        std::cout << "= " << sumper << " <= " << maxEnergy << std::endl;
        std::cout << "dim = " << dim << std::endl;
    }
    eigens.resize(dim);
}


}
