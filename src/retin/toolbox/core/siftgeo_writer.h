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
 * \file siftgeo_writer.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __siftgeo_writer_h__
#define __siftgeo_writer_h__

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <stdexcept>

namespace retin {

template<class T>
class siftgeo_writer {
public:   
    
    static bool write_file(const std::string& file_name,const float* coords,const T* descriptors,size_t feature_dim,size_t feature_count);
    
    siftgeo_writer(std::ostream& out);
    virtual ~siftgeo_writer();    
    
    void        write_feature(const float* coords,const T* descriptors,size_t featureDim,size_t featureCount=1);
    
protected:
    std::ostream& out;
};


inline bool    write_siftgeo_file (const std::string& file_name,const float* coords,const float* descriptors,size_t feature_dim,size_t feature_count) {
    uint8_t* temp = new uint8_t[feature_dim*feature_count];
    for (size_t i=0;i<feature_count*feature_dim;i++) {
        float x = descriptors[i];
        if (x < 0) x = 0;
        else if (x > 255) x = 255;
        temp[i] = uint8_t(x);
    }
    bool res = siftgeo_writer<uint8_t>::write_file(file_name,coords,temp,feature_dim,feature_count);
    delete[] temp;
    return res;
}



template<class T>
bool siftgeo_writer<T>::write_file(const std::string& file_name,const float* coords,const T* descriptors, size_t feature_dim, size_t feature_count) {
    std::ofstream file;
    file.open(file_name.c_str(),std::ios::binary);
    if (file.fail())
        return false;
    try {
        siftgeo_writer<T> writer(file);
        writer.write_feature(coords,descriptors,feature_dim,feature_count);
        file.close();
    }
    catch(...) {
        return false;
    }
    return true;
}

template<class T>
siftgeo_writer<T>::siftgeo_writer(std::ostream& out) : out(out) {
    
}

template<class T>
siftgeo_writer<T>::~siftgeo_writer() {
}

template<class T>
void siftgeo_writer<T>::write_feature(const float* coords,const T* descriptors, size_t featureDim, size_t featureCount) {
    for (size_t i=0;i<featureCount;i++) {
        out.write((const char*)(coords+i*9),9*sizeof(float));
        uint32_t dim = featureDim;
        out.write((const char*)&dim,sizeof(dim));
        if (out.bad())
            throw std::runtime_error("Output stream error");
        out.write((const char*)(descriptors+i*dim),dim*sizeof(T));
        if (out.bad())
            throw std::runtime_error("Output stream error");
    }
}


}

#endif
