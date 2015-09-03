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
 * \file fvec_reader.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __fvec_reader_h__
#define __fvec_reader_h__

#include "string_retin.h"
#include "siftgeo_reader.h"

#include <stdint.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

template<class T>
class fvec_reader {
public:    
    
    static   bool  read_file(const std::string& file_name,T*& buffer,size_t& feature_dim,size_t& feature_count);
    
    fvec_reader(std::istream& in);
    virtual ~fvec_reader();
    
    size_t      get_feature_count();
    size_t      get_feature_dim();
    
    void        read_feature (T* buffer,size_t truncateDim=0);
    void        read_all (T*& buffer,size_t truncateDim=0);
        
protected:
        std::istream& in;
        size_t      featureCount,featureDim;
};

template<class T>
bool  fvec_reader<T>::read_file(const std::string& file_name,T*& buffer,size_t& feature_dim,size_t& feature_count) {
    std::ifstream file;
    file.open(file_name.c_str(),std::ios::binary);
    if (file.fail())
        return false;
    try {
        fvec_reader<T> reader(file);
        reader.read_all(buffer);
        feature_count = reader.get_feature_count();
        feature_dim = reader.get_feature_dim();
        file.close();
    }
    catch(std::exception& ex) {
		std::cerr << ex.what() << std::endl;
        return false;
    }
    return true;
}


template<class T>
fvec_reader<T>::fvec_reader(std::istream& in) : in(in),featureCount(0),featureDim(0) {
    
}

template<class T>
fvec_reader<T>::~fvec_reader() {
}

template<class T>
size_t fvec_reader<T>::get_feature_count() {
    if (featureCount != 0)
        return featureCount;
    size_t dim = get_feature_dim();
    in.seekg(0,std::ios::end);
    if (in.fail())
        throw std::runtime_error("Input stream error 1");
    size_t fileSize = in.tellg();
    size_t rowSize = sizeof(uint32_t) + dim * sizeof(T);
    if ( (fileSize%rowSize) != 0)
        throw std::runtime_error("Invalid file size");
    featureCount = fileSize / rowSize;
    in.seekg(0,std::ios::beg);
    if (in.fail())
        throw std::runtime_error("Input stream error 2");
    return featureCount;
}

template<class T>
size_t fvec_reader<T>::get_feature_dim() {
    if (featureDim != 0)
        return featureDim;
    uint32_t dim = 0;
    in.read((char*)&dim,sizeof(dim));    
    if (in.fail())
        throw std::runtime_error("Input stream error 3");
    in.seekg(0,std::ios::beg);
    if (in.fail())
        throw std::runtime_error("Input stream error 4");
    featureDim = dim;
    return featureDim;
}

template<class T>
void fvec_reader<T>::read_feature (T* buffer,size_t truncateDim) {
    uint32_t dim = 0;
    in.read((char*)&dim,sizeof(dim));    
    if (in.fail())
        throw std::runtime_error("Input stream error 5");
    if (featureDim == 0)
        featureDim = dim;
    else if (featureDim != dim)
        throw std::runtime_error("Invalid feature dim");
    if (truncateDim == 0)
        in.read((char*)buffer,dim*sizeof(T));
    else if (truncateDim > dim)
        throw std::runtime_error("Invalid truncateDim dim");
    else {
        in.read((char*)buffer,truncateDim*sizeof(T));
            if (truncateDim < dim)
                    in.seekg( (dim-truncateDim)*sizeof(T),std::ios::cur );
    }
    if (in.fail())
        throw std::runtime_error("Input stream error 6");
}

template<class T>
void fvec_reader<T>::read_all (T*& buffer,size_t truncateDim) {
    if (buffer)
        throw std::runtime_error("Output buffer in not empty");
    buffer = new T[ get_feature_dim()*get_feature_count() ];
    if (!buffer)
        throw std::runtime_error("Allocation error");
    for (size_t i=0;i<featureCount;i++) 
        read_feature( buffer+i*featureDim,truncateDim );
}

template<class T>
bool	read_bin_desc_file(const std::string& file_name,T*& descriptors,size_t& descriptorDim,size_t& descriptorsCount)
{
    if (descriptors)
        throw std::runtime_error("Output buffer in not empty");
        
	std::ifstream file;
    file.open(file_name.c_str(),std::ios::binary);
    if (file.fail())
        return false;
    try {
		T x = 0;

		file.read((char*)&x,sizeof(x));    
		if (file.fail())
		    throw std::runtime_error("Input stream error 1");
		descriptorsCount = ((x&0xFF)<<8) + (x>>8);

		file.read((char*)&x,sizeof(x));    
		if (file.fail())
		    throw std::runtime_error("Input stream error 2");
		descriptorDim = ((x&0xFF)<<8) + (x>>8);

		descriptors = new T[ descriptorsCount*descriptorDim ];
		file.read((char*)descriptors,descriptorsCount*descriptorDim*sizeof(T));
		if (file.fail() || size_t(file.gcount()) != descriptorsCount*descriptorDim*sizeof(T))
		    throw std::runtime_error("Input stream error 3");
	    file.close();
        
        for (size_t i=0;i<descriptorsCount*descriptorDim;i++) {
            x = descriptors[i];
            descriptors[i] = ((x&0xFF)<<8) + (x>>8);
        }
	}
    catch(...) {
        return false;
    }
    return true;
}

inline bool	read_peter8_desc_file(const std::string& file_name,float*& descriptors,size_t& descriptorDim,size_t& descriptorsCount)
{
    if (descriptors)
        throw std::runtime_error("Output buffer in not empty");
        
	std::ifstream file;
    file.open(file_name.c_str(),std::ios::binary);
    if (file.fail())
        return false;

	file.seekg(0,std::ios::end);
	if (file.fail())
	    throw std::runtime_error("Input stream error 0");
	size_t file_size = file.tellg();
	file.seekg(0,std::ios::beg);

	uint16_t x = 0;

	file.read((char*)&x,sizeof(x));    
	if (file.fail())
	    throw std::runtime_error("Input stream error 1");
	descriptorsCount = x;

	file.read((char*)&x,sizeof(x));    
	if (file.fail())
	    throw std::runtime_error("Input stream error 2");
	descriptorDim = x;

	
	while(1) {
		size_t expected_size = 4+descriptorsCount*(2+descriptorDim);
		if (expected_size == file_size)
			break;
		if (expected_size > file_size) {
		    throw std::runtime_error("Bad file size");
		}
		descriptorsCount += 65536;
	}


	uint8_t* temp = new uint8_t[ descriptorDim+2 ];
	descriptors = new float[ descriptorsCount*descriptorDim ];
	for (size_t i=0;i<descriptorsCount;i++) {
		file.read((char*)temp,descriptorDim+2);
		if (file.fail() || size_t(file.gcount()) != descriptorDim+2)
			throw std::runtime_error("Input stream error 3");
		descriptors[0+i*descriptorDim] = ((uint16_t*)temp)[0];
		descriptors[1+i*descriptorDim] = ((uint16_t*)temp)[1];
		size_t k;
		for (k=2;k<descriptorDim;k++) {
			descriptors[k+i*descriptorDim] = temp[2+k];
		}
	}
	delete[] temp;

	std::streampos pos = file.tellg();
	file.seekg(0,std::ios::end);
	if (file.tellg() != pos) {
		std::cout << "Warning ! Bad file size '" << file_name << "' : " << pos << " != " << file.tellg() << std::endl;
		throw std::runtime_error("Input stream error 4");       
	} 

    file.close(); 
    return true;
}

template<class T>
bool	read_m_desc_file(const std::string& file_name,T*& descriptors,size_t& descriptorDim,size_t& descriptorsCount,char sep=',')
{
    try {
        if (descriptors)
            throw std::runtime_error("Output buffer in not empty");
        
        std::ifstream file;
        char line[0x100000];
        char value[0x100];        
        // count lines;
        file.open(file_name.c_str());
        if (file.fail())
            return false;
        descriptorDim = 0;
        descriptorsCount = 0;
        while (file.good()) {
            line[0] = 0;
            file.getline(line,0x100000);
            if (line[0] == 0 || file.eof())
                break;
            if (file.fail())
                throw std::runtime_error("Input stream error 1");
            if (descriptorsCount == 0) {
                std::istringstream iss(line);
                while(iss.getline(value,0x100,sep))
                    descriptorsCount ++;
            }
            descriptorDim ++;
        }
        file.close();

        if (descriptorDim == 0)
           throw std::runtime_error("Descriptors dimension is zero");
        if (descriptorsCount == 0)
           throw std::runtime_error("No descriptors found");

        // alloc
        descriptors = new T [descriptorDim*descriptorsCount];
        // read data
        size_t j = 0;
        file.open(file_name.c_str());
        while (file.good()) {
            line[0] = 0;
            file.getline(line,0x100000);
            if (file.eof() || line[0] == 0)
                break;
            if (!file.good())
                throw std::runtime_error("Input stream error 2");
            std::istringstream iss(line);
            size_t i = 0;
            while(iss.getline(value,0x100,sep)) {
                double x = strtod(value,NULL);
                descriptors[j+i*descriptorDim] = (T)x;
                i ++;
            }
            j ++;
        }
		file.close();
	}
    catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
    return true;
}

template<class T>
bool	read_text_desc_file(const std::string& file_name,T*& descriptors,size_t& descriptorDim,size_t& descriptorsCount,char sep=',')
{
    try {
        if (descriptors)
            throw std::runtime_error("Output buffer in not empty");

        std::ifstream file;
        char line[0x10000];
        char value[0x100];        
        // count lines;
        file.open(file_name.c_str());
        if (file.fail())
            return false;
        descriptorDim = 0;
        descriptorsCount = 0;
        while (file.good()) {
            line[0] = 0;
            file.getline(line,0x10000);
            if (line[0] == 0)
                break;
            if (file.fail())
                throw std::runtime_error("Input stream error 1");
            if (descriptorDim == 0) {
                std::istringstream iss(line);
                while(iss.getline(value,0x100,sep))
                    descriptorDim ++;
            }
            descriptorsCount ++;
        }
        file.close();

        if (descriptorDim == 0)
           throw std::runtime_error("Descriptors dimension is zero");
        if (descriptorsCount == 0)
           throw std::runtime_error("No descriptors found");

        // alloc
        descriptors = new T [descriptorDim*descriptorsCount];
        // read data
        size_t j = 0;
        file.open(file_name.c_str());
        while (file.good()) {
            line[0] = 0;
            file.getline(line,0x10000);
            if (file.eof() || line[0] == 0)
                break;
            if (!file.good())
                throw std::runtime_error("Input stream error 2");
            std::istringstream iss(line);
            size_t i = 0;
            while(iss.getline(value,0x100,sep)) {
                double x = strtod(value,NULL);
                descriptors[i+j*descriptorDim] = (T)x;
                i ++;
            }
            j ++;
        }
		file.close();
	}
    catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
    return true;
}

inline bool    read_fvec_file (const std::string& file_name,float*& buffer,size_t& feature_dim,size_t& feature_count) {
    if (retin::string_has_file_ext(file_name,".fvec")
     || retin::string_has_file_ext(file_name,".fvecs")) {
        return fvec_reader<float>::read_file(file_name,buffer,feature_dim,feature_count);
    }
    else if (retin::string_has_file_ext(file_name,".hvec8")) {
        uint8_t* temp = NULL;
        if (!fvec_reader<uint8_t>::read_file(file_name,temp,feature_dim,feature_count))
            return false;
        size_t n = feature_dim*feature_count;
        buffer = new float[n];
        for (size_t i=0;i<n;i++) buffer[i] = (float)temp[i];
        delete[] temp;
        return true;
    }
    else if (retin::string_has_file_ext(file_name,".ivec") 
          || retin::string_has_file_ext(file_name,".ivecs")) {
        uint32_t* temp = NULL;
        if (!fvec_reader<uint32_t>::read_file(file_name,temp,feature_dim,feature_count))
            return false;
        size_t n = feature_dim*feature_count;
        buffer = new float[n];
        for (size_t i=0;i<n;i++) buffer[i] = (float)temp[i];
        delete[] temp;
        return true;
    }
    else if (retin::string_has_file_ext(file_name,".siftgeo")) {
        float* coords = NULL;
        bool res = read_siftgeo_file(file_name,coords,buffer,feature_dim,feature_count);
        delete[] coords;
        return res;
    }
    else if (retin::string_has_file_ext(file_name,".m")) {
        return read_m_desc_file<float>(file_name,buffer,feature_dim,feature_count,',');
    }
    else if (retin::string_has_file_ext(file_name,".csv")) {
        return read_text_desc_file<float>(file_name,buffer,feature_dim,feature_count,',');
    }
    else if (retin::string_has_file_ext(file_name,".txt")) {
        return read_text_desc_file<float>(file_name,buffer,feature_dim,feature_count,' ');
    }
    else if (retin::string_has_file_ext(file_name,".bin16")) {
        uint16_t* temp = NULL;
        if (!read_bin_desc_file<uint16_t>(file_name,temp,feature_dim,feature_count))
            return false;
        size_t n = feature_dim*feature_count;
        buffer = new float[n];
        for (size_t i=0;i<n;i++) buffer[i] = (float)temp[i];
        delete[] temp;
        return true;
    }
    else if (retin::string_has_file_ext(file_name,".peter8")) {
		return read_peter8_desc_file(file_name,buffer,feature_dim,feature_count);
    }
    std::cerr << "Unsupported file format " << file_name << std::endl;
    return false;
}


#endif
