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
#ifndef __fvec_reader_h__
#define __fvec_reader_h__

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>

///////////////////////////////////
// BASE CLASS FOR MATRIX READERS //
///////////////////////////////////

template<class T> class matrix_reader {
public:
	matrix_reader(const char* filename, bool bFormatBinary = false) {
		if (!filename)
			return;

		if (bFormatBinary)
			in.open(filename, std::ios::binary);
		else
			in.open(filename);
		if (in.fail())
			throw "Can not open file";
	}
	virtual ~matrix_reader() {
		if (in.is_open())
			in.close();
	}

	virtual size_t get_height() = 0;
	virtual size_t get_width() = 0;

	virtual void read(T*& buffer, size_t truncateDim = 0) = 0;

protected:
	std::ifstream in;
};

/////////////////////////
// GENERAL FVEC READER //
/////////////////////////

template<class T> class general_fvec_reader: public matrix_reader<T> {
public:
	general_fvec_reader(const char* filename);
	virtual ~general_fvec_reader();

	virtual size_t get_height();
	virtual size_t get_width();

	virtual void read(T*& buffer, size_t truncateDim = 0);

protected:
	void read_height();
	void read_width();
	void read_feature(T* buffer, size_t truncateDim = 0);

protected:
	size_t w, h;
};

////////////////////////
// GENERAL BIN READER //
////////////////////////

template<class T> class general_bin_reader: public matrix_reader<T> {
public:
	general_bin_reader(const char* filename);
	virtual ~general_bin_reader();

	virtual size_t get_height();
	virtual size_t get_width();

	virtual void read(T*& buffer, size_t truncateDim = 0);

protected:
	size_t w, h;
};

/////////////////////// READERS /////////////////////////

/////////////////
// FVEC READER //
/////////////////

typedef general_fvec_reader<float> fvec_reader;
typedef general_fvec_reader<double> fvec_reader_double;

//////////////////
// HVEC8 READER //
//////////////////

class hvec8_reader: public matrix_reader<float> {
public:
	hvec8_reader(const char* filename);
	virtual ~hvec8_reader();

	virtual size_t get_height();
	virtual size_t get_width();

	virtual void read(float*& buffer, size_t truncateDim = 0);

protected:
	general_fvec_reader<uint8_t> reader;
};

////////////////
// BIN READER //
////////////////

class bin_reader: public matrix_reader<float> {
public:
	bin_reader(const char* filename);
	virtual ~bin_reader();

	virtual size_t get_height();
	virtual size_t get_width();

	virtual void read(float*& buffer, size_t truncateDim = 0);

protected:
	general_bin_reader<uint16_t> reader;
};

////////////////
// CSV READER //
////////////////

template<class T> class csv_reader: public matrix_reader<T> {
public:
	csv_reader(const char* filename);
	virtual ~csv_reader();

	virtual size_t get_height();
	virtual size_t get_width();

	virtual void read(T*& buffer, size_t truncateDim = 0);

protected:
	size_t w, h;
};

////////////////
// TXT READER //
////////////////

template<class T> class txt_reader: public matrix_reader<T> {
public:
	txt_reader(const char* filename);
	virtual ~txt_reader();

	virtual size_t get_height();
	virtual size_t get_width();

	virtual void read(T*& buffer, size_t truncateDim = 0);

protected:
	size_t w, h;
};

///////////////////
// PETER8 READER //
///////////////////

class peter8_reader: public matrix_reader<float> {
public:
	peter8_reader(const char* filename);
	virtual ~peter8_reader();

	virtual size_t get_height();
	virtual size_t get_width();

	virtual void read(float*& buffer, size_t truncateDim = 0);

protected:
	size_t w, h;
};

///////////////////////
// IDX3 UBYTE READER //
///////////////////////

class idx3_ubyte_reader: public matrix_reader<float> {
public:
	idx3_ubyte_reader(const char* filename);
	virtual ~idx3_ubyte_reader();

	virtual size_t get_height();
	virtual size_t get_width();

	virtual void read(float*& buffer, size_t truncateDim = 0);

protected:
	size_t w, h;
};

/////////////// TEMPLATE METHODS ///////////////

////////////////////////
// csv_reader methods //
////////////////////////

template<class T> csv_reader<T>::csv_reader(const char* filename) :
		matrix_reader<T>(filename, false), w(0), h(0) {
	char line[0x10000];
	char value[0x100];

	// count lines;
	w = 0;
	h = 0;
	while (!this->in.eof()) {
		line[0] = 0;
		this->in.getline(line, 0x10000);
		if (line[0] == 0)
			break;
		if (this->in.fail())
			throw std::runtime_error("Input stream error 1");
		if (w == 0) {
			std::istringstream iss(line);
			while (iss.getline(value, 0x100, ','))
				w++;
		}
		h++;
	}
	this->in.clear();
	this->in.seekg(0, std::ios::beg);
}

template<class T> csv_reader<T>::~csv_reader() {
}

template<class T> size_t csv_reader<T>::get_height() {
	return h;
}
template<class T> size_t csv_reader<T>::get_width() {
	return w;
}

template<class T> void csv_reader<T>::read(T*& data, size_t truncateDim) {
	char line[0x10000];
	char value[0x100];

	size_t j = 0;
	while (this->in.good()) {
		line[0] = 0;
		this->in.getline(line, 0x10000);
		if (line[0] == 0)
			break;
		if (this->in.fail())
			throw std::runtime_error("Input stream error 1");
		std::istringstream iss(line);
		size_t i = 0;
		while (iss.getline(value, 0x100, ',')) {
			double x = strtod(value, NULL);
			data[j + i * w] = (T) x;
			i++;
		}
		j++;
	}
}

////////////////////////
// txt_reader methods //
////////////////////////

template<class T> txt_reader<T>::txt_reader(const char* filename) :
		matrix_reader<T>(filename, false), w(0), h(0) {

	char line[0x10000];
	char value[0x100];

	// count lines;
	w = 0;
	h = 0;
	while (!this->in.eof()) {
		line[0] = 0;
		this->in.getline(line, 0x10000);
		if (line[0] == 0)
			break;
		if (this->in.fail())
			throw "Input stream error 1";
		if (w == 0) {
			std::istringstream iss(line);
			while (iss.getline(value, 0x100, ' '))
				w++;
		}
		h++;
	}
	this->in.clear();
	this->in.seekg(0, std::ios::beg);
}

template<class T> txt_reader<T>::~txt_reader() {
}

template<class T> size_t txt_reader<T>::get_height() {
	return h;
}
template<class T> size_t txt_reader<T>::get_width() {
	return w;
}

template<class T> void txt_reader<T>::read(T*& data, size_t truncateDim) {
	char line[0x10000];
	char value[0x100];

	size_t i = 0;
	while (this->in.good()) {
		line[0] = 0;
		this->in.getline(line, 0x10000);
		if (line[0] == 0)
			break;
		if (this->in.fail())
			throw "Input stream error 1";
		std::istringstream iss(line);
		size_t j = 0;
		while (iss.getline(value, 0x100, ' ')) {
			double x = strtod(value, NULL);
			data[j + i * w] = (T) x;
			j++;
		}
		i++;
	}
}

/////////////////////////
// fvec_reader methods //
/////////////////////////

template<class T> general_fvec_reader<T>::general_fvec_reader(
		const char* filename) :
		matrix_reader<T>(filename, true) {
	h = w = 0;
	read_width();
	read_height();
}
template<class T> general_fvec_reader<T>::~general_fvec_reader() {
}

template<class T> size_t general_fvec_reader<T>::get_height() {
	if (h == 0)
		read_height();
	return h;
}

template<class T> size_t general_fvec_reader<T>::get_width() {
	if (w == 0)
		read_width();
	return w;
}

template<class T> void general_fvec_reader<T>::read_height() {
	size_t dim = get_width();
	this->in.seekg(0, std::ios::end);
	if (this->in.fail())
		throw std::runtime_error("Input stream error 1");
	size_t fileSize = this->in.tellg();
	size_t rowSize = sizeof(uint32_t) + dim * sizeof(T);
	if ((fileSize % rowSize) != 0)
		throw std::runtime_error("Invalid file size");
	h = fileSize / rowSize;
	this->in.seekg(0, std::ios::beg);
	if (this->in.fail())
		throw std::runtime_error("Input stream error 2");
}

template<class T> void general_fvec_reader<T>::read_width() {
	uint32_t dim = 0;
	this->in.read((char*) &dim, sizeof(dim));
	if (this->in.fail())
		throw std::runtime_error("Input stream error 3");
	this->in.seekg(0, std::ios::beg);
	if (this->in.fail())
		throw std::runtime_error("Input stream error 4");
	w = dim;
}

template<class T> void general_fvec_reader<T>::read_feature(T* buffer,
		size_t truncateDim) {
	uint32_t dim = 0;
	this->in.read((char*) &dim, sizeof(dim));
	if (this->in.fail())
		throw std::runtime_error("Input stream error 5");
	if (w == 0)
		w = dim;
	else if (w != dim)
		throw std::runtime_error("Invalid feature dim");
	if (truncateDim == 0)
		this->in.read((char*) buffer, dim * sizeof(T));
	else if (truncateDim > dim)
		throw std::runtime_error("Invalid truncateDim dim");
	else {
		this->in.read((char*) buffer, truncateDim * sizeof(T));
		if (truncateDim < dim)
			this->in.seekg((dim - truncateDim) * sizeof(T), std::ios::cur);
	}
	if (this->in.fail())
		throw std::runtime_error("Input stream error 6");
}

template<class T> void general_fvec_reader<T>::read(T*& buffer,
		size_t truncateDim) {
	if (!buffer)
		buffer = new T[get_width() * get_height()];
	for (size_t i = 0; i < h; i++)
		read_feature(buffer + i * w, truncateDim);
}

////////////////////////////////
// general_bin_reader methods //
////////////////////////////////

template<class T> general_bin_reader<T>::general_bin_reader(
		const char* filename) :
		matrix_reader<T>(filename, true) {
	h = w = 0;
	T x = 0;

	this->in.read((char*) &x, sizeof(x));
	if (this->in.fail())
		throw std::runtime_error("Input stream error 1");
	h = ((x & 0xFF) << 8) + (x >> 8);

	this->in.read((char*) &x, sizeof(x));
	if (this->in.fail())
		throw std::runtime_error("Input stream error 2");
	w = ((x & 0xFF) << 8) + (x >> 8);
}

template<class T> general_bin_reader<T>::~general_bin_reader() {
}

template<class T> size_t general_bin_reader<T>::get_height() {
	return h;
}
template<class T> size_t general_bin_reader<T>::get_width() {
	return w;
}

template<class T> void general_bin_reader<T>::read(T*& data,
		size_t truncateDim) {
	T x = 0;

	this->in.read((char*) data, h * w * sizeof(T));
	if (this->in.fail() || size_t(this->in.gcount()) != h * w * sizeof(T))
		throw std::runtime_error("Input stream error 3");

	for (size_t i = 0; i < h * w; i++) {
		x = data[i];
		data[i] = ((x & 0xFF) << 8) + (x >> 8);
	}
}
#endif
