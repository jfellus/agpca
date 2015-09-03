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

#include "matrix_writer.h"



///////////////////
// peter8_writer //
///////////////////

void peter8_writer::write(const float* buffer, size_t w, size_t h) {
		uint16_t x = 0;

		// header
		x = h;
		out.write((const char*)&x,sizeof(x));
		if (out.fail())
			throw std::runtime_error("Input stream error 1");

		x = w;
		out.write((const char*)&x,sizeof(x));
		if (out.fail())
			throw std::runtime_error("Input stream error 2");

		// descriptors
		for (size_t i=0;i<h;i++) {
			const float* row = buffer + i*w;
			x = row[0];
			out.write((const char*)&x,sizeof(x));
			x = row[1];
			out.write((const char*)&x,sizeof(x));
			for (size_t k=2;k<w;k++) {
				uint8_t b = row[k];
				out.write((const char*)&b,sizeof(b));
			}
			if (out.fail())
				throw std::runtime_error("Input stream error 3");
		}

}


//////////////////
// hvec8_writer //
//////////////////

void hvec8_writer::write(const float* buffer, size_t w, size_t h) {
	uint8_t* temp = 0;
	size_t n = w*h;
	temp = new uint8_t[n];
	for (size_t i=0;i<n;i++) {
		float x = buffer[i];
		if (x < 0 || x > 255) {
			std::cout << "Index: " << i << ", value: " << x << std::endl;
			throw std::runtime_error("Values outside [0,255]");
		}
		temp[i] = (uint8_t)x;
	}
	writer.write(temp, w, h);
	delete[] temp;
}


//////////////////
// ivecs_writer //
//////////////////

void ivecs_writer::write(const float* buffer, size_t w, size_t h) {
    uint32_t* temp = NULL;
    size_t n = w*h;
    temp = new uint32_t[n];
    for (size_t i=0;i<n;i++) {
        float x = buffer[i];
        if (x < 0 || x > 0xFFFFFFFF) {
            std::cout << "Index: " << i << ", value: " << x << std::endl;
            throw std::runtime_error("Values outside range");
        }
        temp[i] = (uint32_t)x;
    }
    writer.write(temp, w, h);
    delete[] temp;
}

