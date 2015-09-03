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

#ifndef MATRIX_WRITER_H_
#define MATRIX_WRITER_H_

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>

///////////////////////////////////
// BASE CLASS FOR MATRIX WRITERS //
///////////////////////////////////

template<class T> class matrix_writer {
public:
	matrix_writer(const char* filename, bool bFormatBinary = false) {
		if (!filename)
			return;

		if (bFormatBinary)
			out.open(filename, std::ios::binary);
		else
			out.open(filename);
		if (out.fail())
			throw "Can not open file for writing";
	}
	virtual ~matrix_writer() {
		if (out.is_open())
			out.close();
	}

	virtual void write(const T* buffer, size_t w, size_t h) = 0;

protected:
	std::ofstream out;
};

/////////////////////////
// GENERAL FVEC WRITER //
/////////////////////////

template<class T> class general_fvec_writer: public matrix_writer<T> {
public:
	general_fvec_writer(const char* filename) :
			matrix_writer<T>(filename, true) {
	}
	virtual ~general_fvec_writer() {
	}

	virtual void write(const T* buffer, size_t w, size_t h);
};

/////////////////////// WRITERS /////////////////////////

/////////////////
// FVEC WRITER //
/////////////////

typedef general_fvec_writer<float> fvec_writer;

//////////////////
// HVEC8 WRITER //
//////////////////

class hvec8_writer: public matrix_writer<float> {
public:
	hvec8_writer(const char* filename) :
			matrix_writer<float>(0), writer(filename) {
	}
	virtual ~hvec8_writer() {
	}

	virtual void write(const float* buffer, size_t w, size_t h);

protected:
	general_fvec_writer<uint8_t> writer;
};

//////////////////
// IVECS WRITER //
//////////////////

class ivecs_writer: public matrix_writer<float> {
public:
	ivecs_writer(const char* filename) :
			matrix_writer<float>(0), writer(filename) {
	}
	virtual ~ivecs_writer() {
	}

	virtual void write(const float* buffer, size_t w, size_t h);

protected:
	general_fvec_writer<uint32_t> writer;
};

////////////////
// CSV WRITER //
////////////////

template<class T> class csv_writer: public matrix_writer<T> {
public:
	csv_writer(const char* filename) :
			matrix_writer<T>(filename, false) {
	}
	virtual ~csv_writer() {
	}

	virtual void write(const T* buffer, size_t w, size_t h);
};

////////////////
// TXT WRITER //
////////////////

template<class T> class txt_writer: public matrix_writer<T> {
public:
	txt_writer(const char* filename) :
			matrix_writer<T>(filename, false) {
	}
	virtual ~txt_writer() {
	}

	virtual void write(const T* buffer, size_t w, size_t h);
};

///////////////////
// PETER8 WRITER //
///////////////////

class peter8_writer: public matrix_writer<float> {
public:
	peter8_writer(const char* filename) :
			matrix_writer<float>(filename, true) {
	}
	virtual ~peter8_writer() {
	}

	virtual void write(const float* buffer, size_t w, size_t h);
};

///////////////// TEMPLATES METHODS /////////////////

////////////////
// csv_writer //
////////////////

template<class T> void csv_writer<T>::write(const T* buffer, size_t w,
		size_t h) {
	for (size_t k = 0; k < w; k++) {
		this->out << buffer[k];
		for (size_t i = 1; i < h; i++) {
			this->out << "," << buffer[k + i * w];
		}
		this->out << '\n';
	}
}

////////////////
// csv_writer //
////////////////

template<class T> void txt_writer<T>::write(const T* buffer, size_t w,
		size_t h) {
	for (size_t i = 0; i < h; i++) {
		this->out << buffer[i * w];
		for (size_t j = 1; j < w; j++) {
			this->out << " " << buffer[j + i * w];
		}
		this->out << '\n';
	}
}

/////////////////////////
// general_fvec_writer //
/////////////////////////

template<class T> void general_fvec_writer<T>::write(const T* buffer, size_t w,
		size_t h) {
	for (size_t i = 0; i < h; i++) {
		uint32_t dim = w;
		this->out.write((const char*) &dim, sizeof(dim));
		if (this->out.bad())
			throw std::runtime_error("Output stream error");
		this->out.write((const char*) (buffer + i * dim), dim * sizeof(T));
		if (this->out.bad())
			throw std::runtime_error("Output stream error");
	}
}

#endif
