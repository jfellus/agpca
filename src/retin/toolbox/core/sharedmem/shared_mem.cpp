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

#include "shared_mem.h"

void* create_shared_mem(key_t key, unsigned long size, int* _shmid_out) {
	*_shmid_out = shmget(key, size, IPC_CREAT | IPC_EXCL | 0666);
	if(*_shmid_out < 0) {
		std::cerr << "Maybe too large..." << size << "\n";
		throw "Unable to allocate shared memory";
	}

	void* data = shmat(*_shmid_out, 0, 0);
	if(data == (char*)-1) throw "Unable to attach to shared memory";

	return data;
}

void delete_shared_mem(void* ptr, int _shmid, const char* id) {
	detach_shared_mem(ptr);
	unlink(id);
	if(shmctl(_shmid, IPC_RMID, 0) < 0) throw "Unable to delete shared memory";
}

void* attach_shared_mem(key_t key, int* _shmid_out) {
	*_shmid_out = shmget(key, 0, 0666);
	if(*_shmid_out < 0) throw "No shared memory found with this id";

	void* data = shmat(*_shmid_out, 0, 0);
	if(data == (char*)-1) throw "Unable to attach to shared memory";

	return data;
}

void detach_shared_mem(void* ptr) {
	if(shmdt(ptr)<0) throw "Unable to detach shared memory (given pointer does not point to shared segment";
}


