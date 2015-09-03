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

#include "SharedMatrix.h"
#include "sharedmem/shared_mem.h"
#include "sharedmem/matrix_reader.h"
#include "sharedmem/matrix_writer.h"
#include <string.h>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>


using namespace std;
using namespace retin;

static bool __SHARED_MATRIX_DISABLED = false;

///////////
// UTILS //
///////////

static string str_replace(string subject, const string& search, const string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != string::npos) {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
    }
    return subject;
}

static void touch(const char* file) {
	FILE* f = fopen(file,"w"); if(f) fclose(f);
}

static key_t create_key(const string& file) {
	touch(file.c_str());
	return ftok(file.c_str(), 1);
}

static const char* getcwd() {
	static char* c = NULL;
	if(c==NULL) c = new char[512];
	if(!getcwd(c, 512)) {};
	return c;
}

static string realpath(const char* file) {
	string realfile = getcwd();
	if(file[0]=='/') realfile = file;
	else {realfile += "/"; realfile += file;}
	char* rp = realpath(realfile.c_str(),NULL);
	if(rp) {
		realfile = rp;
		free(rp);
	}
	return realfile;
}

static string create_id(const char* file) {
	if(system("mkdir --parents /var/tmp/shared_matrices/")){}
	return string("/var/tmp/shared_matrices/") + str_replace(realpath(file),"/","__");
}

namespace shared_matrices {

////////////
// MATRIX //
////////////

/** Creates a new shared matrix from a file */
void Matrix::load(const char* file, int flags) {
	if(__SHARED_MATRIX_DISABLED) flags = MATRIX_LOCAL;

	if((flags & MATRIX_LOCAL) && (flags & MATRIX_CREATE))
		throw "MATRIX_LOCAL and MATRIX_CREATE are incompatible";

	if((flags & MATRIX_ATTACH) && attach(file)) {
		this->bShared = true;
		return;
	}

	this->bShared = ((flags & MATRIX_LOCAL)==0);

	if(flags & MATRIX_CREATE) key = create_key(id);
	else key = -1;

	if(flags & MATRIX_LOCAL || flags & MATRIX_CREATE) {
		read(file);
		if(!data) throw "Unrecognized Matrix file format";
	} else data = NULL;
}

/** Attach to an already existing shared matrix */
bool Matrix::attach(const char* file) {
	width=height=0;
	data = 0;
	bOwner = false;

	this->id = create_id(file);

	key = ftok(id.c_str(), 1);
	if(key==-1) return false;

	size_t* _data = NULL;
	try{ _data = (size_t*)attach_shared_mem(key, &_shmid); }catch(...){}
	if(!_data) return false;

	width = *_data++;
	height = *_data++;
	data = (float*)_data;
	bShared = true;

	cout << "Found shared matrix : " << file << "\n";
	return true;
}

void Matrix::detach() {
	 detach_shared_mem(((size_t*)data)-2);
}

void Matrix::delete_shared() {
	delete_shared_mem(((size_t*)data)-2, _shmid, this->id.c_str());
}

void Matrix::create_ref(float* data, size_t width, size_t height) {
	this->width = width;
	this->height = height;
	this->data = data;
}

/** Creates a new empty shared matrix of dimensions  width x height */
void Matrix::create(size_t width, size_t height, const char* shared_id) {
	if(__SHARED_MATRIX_DISABLED) shared_id = NULL;

	if(shared_id) {
		bShared = true;
		if(attach(shared_id)) {
			if(this->width==width && this->height==height) 	return;
			cout << "But matrix dimensions does not agree : reshaping ... ok\n";
			delete_shared();
		}
		key = create_key(id);
		try {
			size_t* _data = (size_t*)create_shared_mem(key, width*height*sizeof(float) + 2*sizeof(size_t), &_shmid);
			*_data++ = width;
			*_data++ = height;
			data = (float*)_data;
			//this->clear();
		} catch(...) {
			data = 0;
			cerr << "Couldn't allocate shared matrix " << shared_id << "\n";
			unlink(id.c_str());
			throw "Couldn't allocate shared matrix";
		}
	} else {
		bShared = false;
		data = new float[width*height];
		//this->clear();
	}

	this->width = width; this->height = height;
}

/** if this instance is the owner of the shared content, delete the whole shared matrix
 *  else destroy only the local reference to the shared matrix */

/** if this instance is the owner of the shared content, delete the whole shared matrix
 *  else destroy only the local reference to the shared matrix */
Matrix::~Matrix() {
	dealloc();
}

void Matrix::dealloc() {
	fflush(stdout);
	if(!data) return;
	if(!bShared) {
		fflush(stdout);
		delete data;
	} else {
		if(bOwner) delete_shared();
		else detach();
	}
}

void Matrix::realloc(size_t w, size_t h, const char* shared_id) {
	dealloc();
	create(w,h,shared_id);
}



bool Matrix::read(const char* file) {
	matrix_reader<float>* reader = 0;
	std::string file_name = file;
	if (string_has_file_ext(file_name,".fvec")) reader = new fvec_reader(file);
    else if (string_has_file_ext(file_name,".hvec8")) reader = new hvec8_reader(file);
    else if (string_has_file_ext(file_name,".csv") || string_has_file_ext(file_name,".m")) reader = new csv_reader<float>(file);
    else if (string_has_file_ext(file_name,".txt")) reader = new txt_reader<float>(file);
    else if (string_has_file_ext(file_name,".bin16")) reader = new bin_reader(file);
    else if (string_has_file_ext(file_name,".peter8")) reader = new peter8_reader(file);
    else if (string_has_file_ext(file_name,".idx3-ubyte")) reader = new idx3_ubyte_reader(file);

	if(!reader) throw "No reader for this file format";

	height = reader->get_height();
	width = reader->get_width();

	if(bShared) {
		size_t* _data = (size_t*)create_shared_mem(key, width*height*sizeof(float) + 2*sizeof(size_t), &_shmid);
		*_data++ = width;
		*_data++ = height;
		data = (float*)_data;
	} else {
		data = new float[width*height];
	}

	fflush(stdout);

    reader->read(data);

	delete reader;

	return true;
}

bool Matrix::write(const char* file) const {
	matrix_writer<float>* writer = 0;
	std::string file_name = file;
	if (string_has_file_ext(file_name,".fvec")) writer = new fvec_writer(file);
	else if (string_has_file_ext(file_name,".hvec8")) writer = new hvec8_writer(file);
	else if (string_has_file_ext(file_name,".csv") || string_has_file_ext(file_name,".m")) writer = new csv_writer<float>(file);
	else if (string_has_file_ext(file_name,".txt")) writer = new txt_writer<float>(file);
	else if (string_has_file_ext(file_name,".ivecs")) writer = new ivecs_writer(file);
	else if (string_has_file_ext(file_name,".peter8")) writer = new peter8_writer(file);

	if(!writer) throw "No writer for this file format";

	writer->write(data, width, height);
	delete writer;
	return true;
}


void Matrix::dump(size_t nbrows, size_t nbcols) {
	if(nbrows==0) nbrows = height;
	if(nbcols==0) nbcols = width;

	for(size_t i=0; i<nbrows; i++) {
		if(i) putc('\n',stdout);
		for(size_t j=0; j<nbcols; j++) {
			if(j) putc(' ',stdout);
			printf("%f", data[i*width+j]);
		}
		fflush(stdout);
	}
	putc('\n',stdout);
}

bool Matrix::free_shared(const char* file) {
	bool b = false;
	Matrix* m = new Matrix(file, MATRIX_ATTACH);
	if(*m) {b=true; cout << "Free shared matrix : " << file << endl;}
	m->set_owner();
	delete m;
	return b;
}

void Matrix::clear() {
	if(!data) return;
	memset(data, 0, width*height*sizeof(float));
}

void Matrix::operator=(const Matrix& m) {memcpy(this->data, m.data, sizeof(float)*width*height);}


}

void SHARED_MATRIX_DISABLE() {
	__SHARED_MATRIX_DISABLED = true;
}
