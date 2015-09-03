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

#ifndef SHAREDMATRIXDOUBLE_H_
#define SHAREDMATRIXDOUBLE_H_

#include <stdlib.h>
#include <string>
#include <string.h>
#include "string_retin.h"

using namespace std;

// MATRIX LOADING FLAGS
#define MATRIX_CREATE 0x0001
#define MATRIX_LOCAL 0x0010
#define MATRIX_ATTACH 0x0100

namespace shared_matrices {

class MatrixDouble {
public:
	size_t width, height;
private:
	bool bOwner, bShared;
	std::string id;

	int _shmid;
	key_t key;
public:
	double* data;

public:
	MatrixDouble() {width=height=0; data=0; bOwner=bShared=false; _shmid=-1;}
	MatrixDouble(const string& file, int flags = MATRIX_LOCAL | MATRIX_ATTACH) {load(file.c_str(), flags);}
	MatrixDouble(const char* file, int flags = MATRIX_LOCAL | MATRIX_ATTACH) {load(file, flags);}
	MatrixDouble(size_t w, size_t h, const char* shared_id = NULL) {create(w,h,shared_id);}
	MatrixDouble(size_t w, size_t h, const string& shared_id) {create(w,h,shared_id.c_str());}

	/** if this instance is the owner of the shared content, delete the whole shared matrix
	 *  else destroy only the local reference to the shared matrix */
	virtual ~MatrixDouble();


	void dealloc();
	void realloc(size_t w, size_t h, const char* shared_id = NULL);



	size_t get_width() const {return width;}
	size_t get_height() const {return height;}

	bool is_allocated() const {return data!=0;}
	bool is_owner() const {return bOwner;}
	bool is_shared() const {return bShared;}


	double* get_row(size_t i) {return data + i*width;}

	/** Make this instance the owner of the shared matrix (i.e. when destroyed, the shared matrix itself
	 * is deleted */
	void set_owner(bool bOwner = true) {
		this->bOwner = bOwner;
	}

	bool read(const char* file);
	bool write(const char* file) const;
	inline bool save(const char* file) const {return this->write(file);}
	inline bool save(const string& file) const {return this->write(file.c_str());}

	void dump(size_t nbrows=0, size_t nbcols=0);

	void clear();

	operator bool() const { return is_allocated(); }
	operator double*() const {return data; }
	double& operator[](int i) {return data[i]; }
	void operator=(const MatrixDouble& m);
	void operator=(const int* data) {for(int i=0; i<width*height; i++) this->data[i] = data[i];}
	void operator=(double x) {for(int i=0; i<width*height; i++) this->data[i] = x;}


	static bool free_shared(const char* file);
	static bool free_shared(const string& file) {return free_shared(file.c_str());}

	void create(size_t w, size_t h, const char* shared_id = NULL);
	void create_ref(double* data, size_t w, size_t h);
	void load(const char* file, int flags = MATRIX_LOCAL);

protected:
	bool attach(const char* file);
	void detach();
	void delete_shared();
};

}

#endif /* SHAREDMATRIXDOUBLE_H_ */
