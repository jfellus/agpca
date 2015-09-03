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
 * \file siftgeo_reader.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __siftgeo_reader_h__
#define __siftgeo_reader_h__

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <stdexcept>

template<class T>
class siftgeo_reader {
public:    
    
    static   bool  read_file(const std::string& file_name,float*& coords,T*& descriptors,size_t& feature_dim,size_t& feature_count);
    
    siftgeo_reader(std::istream& in);
    virtual ~siftgeo_reader();
    
    size_t      get_feature_count();
    size_t      get_feature_dim();
    
    void        read_feature (float* coords,T* descriptors,size_t truncateDim=0);
    void        read_all (float*& coords,T*& descriptors,size_t truncateDim=0);
        
protected:
        std::istream& in;
        size_t      featureCount,featureDim;
};

inline bool    read_siftgeo_file (const std::string& file_name,float*& coords,float*& descriptors,size_t& feature_dim,size_t& feature_count) {
    std::ifstream in;
    in.open(file_name.c_str(),std::ios::binary);
    if (!in)
        return false;
    in.seekg(0,std::ios::end);
    size_t featureDim = 128;
    size_t fileSize = in.tellg();
    size_t rowSize = 9 * 4 + 1 * 4 + featureDim * sizeof(uint8_t);
    if ( (fileSize%rowSize) != 0)
        return false;
    size_t featureCount = fileSize / rowSize;
    in.seekg(0,std::ios::beg);
    if (!in)
        return false;
    uint8_t* temp = new uint8_t[fileSize];
    in.read((char*)temp,fileSize);
    if (!in) {
        return false;
    }
    in.close();
    
    feature_dim = featureDim;
    feature_count = featureCount;
    coords = new float[9*featureCount];
    descriptors = new float[featureDim*featureCount];
    for (size_t i=0;i<featureCount;i++) {
        float* ptr = (float*)(temp+i*rowSize);
        for (size_t k=0;k<9;k++) {
            coords[k+i*9] = ptr[k];
        }
        if ( ((uint32_t*)ptr)[9] != featureDim )
            return false;
        for (size_t k=0;k<featureDim;k++)
            descriptors[k+i*featureDim] = temp[10*4+k+i*rowSize];
    }    

    delete[] temp;
    return true;
}



template<class T>
bool  siftgeo_reader<T>::read_file(const std::string& file_name,float*& coords,T*& descriptors,size_t& feature_dim,size_t& feature_count) {
    std::ifstream file;
    file.open(file_name.c_str(),std::ios::binary);
    if (file.fail())
        return false;
    try {
        siftgeo_reader<T> reader(file);
        reader.read_all(coords,descriptors);
        feature_count = reader.get_feature_count();
        feature_dim = reader.get_feature_dim();
        file.close();
    }
    catch(...) {
        return false;
    }
    return true;
}


template<class T>
siftgeo_reader<T>::siftgeo_reader(std::istream& in) : in(in),featureCount(0),featureDim(128) {
    
}

template<class T>
siftgeo_reader<T>::~siftgeo_reader() {
}

template<class T>
size_t siftgeo_reader<T>::get_feature_count() {
    if (featureCount != 0)
        return featureCount;
    in.seekg(0,std::ios::end);
    if (in.fail())
        throw std::runtime_error("Input stream error");
    size_t fileSize = in.tellg();
    size_t rowSize = 9 * 4 + 1 * 4 + 128 * sizeof(T);
    if ( (fileSize%rowSize) != 0)
        throw std::runtime_error("Invalid file size");
    featureCount = fileSize / rowSize;
    in.seekg(0,std::ios::beg);
    if (in.fail())
        throw std::runtime_error("Input stream error");
    return featureCount;
}

template<class T>
size_t siftgeo_reader<T>::get_feature_dim() {
    return featureDim;
}

template<class T>
void siftgeo_reader<T>::read_feature (float* coords,T* descriptors,size_t truncateDim) {
    in.read((char*)coords,9*sizeof(float));
    if (in.fail())
        throw std::runtime_error("Input stream error");
    uint32_t dim = 0;
    in.read((char*)&dim,sizeof(dim));    
    if (in.fail())
        throw std::runtime_error("Input stream error");
    if (featureDim != dim) 
        throw std::runtime_error("Invalid feature dim");
    if (truncateDim == 0)
        in.read((char*)descriptors,dim*sizeof(T));
    else if (truncateDim > dim)
            throw std::runtime_error("Invalid truncateDim dim");
    else {
        in.read((char*)descriptors,truncateDim*sizeof(T));
        if (truncateDim < dim)
            in.seekg( (dim-truncateDim)*sizeof(T),std::ios::cur );
    }
    if (in.fail())
        throw std::runtime_error("Input stream error");
}

template<class T>
void siftgeo_reader<T>::read_all (float*& coords,T*& descriptors,size_t truncateDim) {
    if (descriptors)
        throw std::runtime_error("Output descriptors in not empty");
    size_t feature_dim = get_feature_dim();
    size_t feature_count = get_feature_count();
    coords = new float[9*feature_count];
    if (!coords)
        throw std::runtime_error("Allocation error");
    descriptors = new T[feature_dim*feature_count];
    if (!descriptors)
        throw std::runtime_error("Allocation error");
    for (size_t i=0;i<featureCount;i++) 
        read_feature( coords+i*9,descriptors+i*featureDim,truncateDim );
}

#endif
