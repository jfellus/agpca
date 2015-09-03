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



#include "matrix_reader.h"


static uint32_t change_endian(uint32_t x) {
	unsigned char* data = (unsigned char*)&x;
	return ((uint32_t)data[0]<<24) | ((uint32_t)data[1]<<16) | ((uint32_t)data[2]<<8) | ((uint32_t)data[3]);
}


///////////////////////////
// peter8_reader methods //
///////////////////////////

peter8_reader::peter8_reader(const char* filename) : matrix_reader<float>(filename, true),w(0),h(0) {
	this->in.seekg(0,std::ios::end);
	if (this->in.fail())
	    throw std::runtime_error("Input stream error 0");
	size_t file_size = this->in.tellg();
	this->in.seekg(0,std::ios::beg);

	uint16_t x = 0;

	this->in.read((char*)&x,sizeof(x));
	if (this->in.fail())
	    throw std::runtime_error("Input stream error 1");
	h = x;

	this->in.read((char*)&x,sizeof(x));
	if (this->in.fail())
	    throw std::runtime_error("Input stream error 2");
	w = x;


	while(1) {
		size_t expected_size = 4+h*(2+w);
		if (expected_size == file_size)
			break;
		if (expected_size > file_size) {
		    throw std::runtime_error("Bad file size");
		}
		h += 65536;
	}

}

peter8_reader::~peter8_reader() {}

size_t peter8_reader::get_height() {return h;}
size_t peter8_reader::get_width() {return w;}

void peter8_reader::read (float*& data,size_t truncateDim) {
	uint8_t* temp = new uint8_t[ w+2 ];
	for (size_t i=0;i<h;i++) {
		this->in.read((char*)temp,w+2);
		if (this->in.fail() || size_t(this->in.gcount()) != w+2)
			throw std::runtime_error("Input stream error 3");
		data[0+i*w] = ((uint16_t*)temp)[0];
		data[1+i*w] = ((uint16_t*)temp)[1];
		size_t k;
		for (k=2;k<w;k++) {
			data[k+i*w] = temp[2+k];
		}
	}
	delete[] temp;

	std::streampos pos = this->in.tellg();
	this->in.seekg(0,std::ios::end);
	if (this->in.tellg() != pos) {
		std::cout << "Warning ! Bad file size " << pos << " != " << this->in.tellg() << std::endl;
		throw std::runtime_error("Input stream error 4");
	}
}

///////////////////////////////
// IDX3 UBYTE READER METHODS //
///////////////////////////////

idx3_ubyte_reader::idx3_ubyte_reader(const char* filename) : matrix_reader<float>(filename, true),w(0),h(0) {
	if (this->in.fail())    throw std::runtime_error("Input stream error 0");

	uint32_t x = 0;

	this->in.read((char*)&x,sizeof(x));
	if (this->in.fail())    throw std::runtime_error("Input stream error 1");

	this->in.read((char*)&x,sizeof(x));
	h = change_endian(x);

	this->in.read((char*)&x,sizeof(x));
	w = change_endian(x);
	this->in.read((char*)&x,sizeof(x));
	w *= change_endian(x);

	std::cout << "D=" << w << "\n";
}

idx3_ubyte_reader::~idx3_ubyte_reader() {}

size_t idx3_ubyte_reader::get_height() {return h;}
size_t idx3_ubyte_reader::get_width() {return w;}

void idx3_ubyte_reader::read (float*& data,size_t truncateDim) {
	for(int i=0; i<w*h; i++) {
		uint8_t temp;
		in.read((char*)&temp, sizeof(temp));
		data[i] = temp/255.0f;
	}
}


//////////////////////////
// hvec8_reader methods //
//////////////////////////


hvec8_reader::hvec8_reader(const char* filename) : matrix_reader<float>(0),reader(filename) {}
hvec8_reader::~hvec8_reader() {}

size_t hvec8_reader::get_height() {return reader.get_height();}
size_t hvec8_reader::get_width() {return reader.get_width();}

void hvec8_reader::read (float*& data,size_t truncateDim) {
	uint8_t* temp = NULL;
	reader.read(temp, truncateDim);
	size_t n = get_width()*get_height();
	for (size_t i=0;i<n;i++) data[i] = (float)temp[i];
	delete[] temp;
}

////////////////////////
// bin_reader methods //
////////////////////////


bin_reader::bin_reader(const char* filename) : matrix_reader<float>(0),reader(filename) {}
bin_reader::~bin_reader() {}

size_t bin_reader::get_height() {return reader.get_height();}
size_t bin_reader::get_width() {return reader.get_width();}

void bin_reader::read (float*& data,size_t truncateDim) {
	uint16_t* temp = NULL;
	reader.read(temp, truncateDim);
	size_t n = get_width()*get_height();
	for (size_t i=0;i<n;i++) data[i] = (float)temp[i];
	delete[] temp;
}





