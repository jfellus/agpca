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
 * \file fvec_writer.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __fvec_writer_h__
#define __fvec_writer_h__

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <stdexcept>

namespace retin {

template<class T>
class fvec_writer {
public:   
    
    static bool write_file(const std::string& file_name,const T* buffer,size_t feature_dim,size_t feature_count);
    
    fvec_writer(std::ostream& out);
    virtual ~fvec_writer();    
    
    void        write_feature(const T* buffer,size_t featureDim,size_t featureCount=1);
    
protected:
    std::ostream& out;
};


template<class T>
bool fvec_writer<T>::write_file(const std::string& file_name,const T* buffer, size_t feature_dim, size_t feature_count) {
    std::ofstream file;
    file.open(file_name.c_str(),std::ios::binary);
    if (file.fail())
        return false;
    try {
        fvec_writer<T> writer(file);
        writer.write_feature(buffer,feature_dim,feature_count);
        file.close();
    }
    catch(...) {
        return false;
    }
    return true;
}

template<class T>
fvec_writer<T>::fvec_writer(std::ostream& out) : out(out) {
    
}

template<class T>
fvec_writer<T>::~fvec_writer() {
}

template<class T>
void fvec_writer<T>::write_feature(const T* buffer, size_t featureDim, size_t featureCount) {
    for (size_t i=0;i<featureCount;i++) {
        uint32_t dim = featureDim;
        out.write((const char*)&dim,sizeof(dim));
        if (out.bad())
            throw std::runtime_error("Output stream error");
        out.write((const char*)(buffer+i*dim),dim*sizeof(T));
        if (out.bad())
            throw std::runtime_error("Output stream error");
    }
}

inline bool write_peter8_desc_file(const std::string& file_name,const float* buffer, size_t feature_dim, size_t feature_count) {
    std::ofstream file;
    file.open(file_name.c_str(),std::ios::binary);
    if (file.fail())
        return false;
    try {
		uint16_t x = 0;

		// header
		x = feature_count;
		file.write((const char*)&x,sizeof(x));
		if (file.fail())
			throw std::runtime_error("Input stream error 1");

		x = feature_dim;
		file.write((const char*)&x,sizeof(x));
		if (file.fail())
			throw std::runtime_error("Input stream error 2");

		// descriptors
		for (size_t i=0;i<feature_count;i++) {
			const float* row = buffer + i*feature_dim;
			x = row[0];
			file.write((const char*)&x,sizeof(x));
			x = row[1];
			file.write((const char*)&x,sizeof(x));
			for (size_t k=2;k<feature_dim;k++) {
				uint8_t b = row[k];
				file.write((const char*)&b,sizeof(b));
			}
			if (file.fail())
				throw std::runtime_error("Input stream error 3");
		}

        file.close();
    }
    catch(...) {
        return false;
    }
    return true;
}

inline bool write_m_desc_file(const std::string& file_name,const float* buffer, size_t feature_dim, size_t feature_count) {
    std::ofstream file;
    file.open(file_name.c_str());
    if (file.fail())
        return false;
    try {
        for (size_t k=0;k<feature_dim;k++) {
            file << buffer[k];
            for (size_t i=1;i<feature_count;i++) {
                file << "," << buffer[k+i*feature_dim];
            }
            file << std::endl;
        }
        file.close();
    }
    catch(...) {
        return false;
    }
    return true;
}

inline bool    write_fvec_file (const std::string& file_name,float* buffer,size_t feature_dim,size_t feature_count) {
    if (file_name.rfind(".fvec") != std::string::npos) {
        return fvec_writer<float>::write_file(file_name,buffer,feature_dim,feature_count);
    }
    else if (file_name.rfind(".hvec8") != std::string::npos) {
        uint8_t* temp = NULL;
        size_t n = feature_dim*feature_count;
        temp = new uint8_t[n];
        for (size_t i=0;i<n;i++) {
            float x = buffer[i];
            if (x < 0 || x > 255) {
                std::cout << "Index: " << i << ", value: " << x << std::endl;
                throw std::runtime_error("Values outside [0,255]");
            }
            temp[i] = (uint8_t)x;
        }
        if (!fvec_writer<uint8_t>::write_file(file_name,temp,feature_dim,feature_count))
            return false;
        delete[] temp;
        return true;
    }
    else if (file_name.rfind(".ivecs") != std::string::npos) {
        uint32_t* temp = NULL;
        size_t n = feature_dim*feature_count;
        temp = new uint32_t[n];
        for (size_t i=0;i<n;i++) {
            float x = buffer[i];
            if (x < 0 || x > 0xFFFFFFFF) {
                std::cout << "Index: " << i << ", value: " << x << std::endl;
                throw std::runtime_error("Values outside range");
            }
            temp[i] = (uint32_t)x;
        }
        if (!fvec_writer<uint32_t>::write_file(file_name,temp,feature_dim,feature_count))
            return false;
        delete[] temp;
        return true;
    }
	else if (file_name.rfind(".peter8") != std::string::npos) {
		return write_peter8_desc_file(file_name,buffer,feature_dim,feature_count);
	}
	else if (file_name.rfind(".m") != std::string::npos) {
		return write_m_desc_file(file_name,buffer,feature_dim,feature_count);
	}
    return false;
}


}

#endif
