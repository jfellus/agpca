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
 * \file auto_array_ptr.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __retin_auto_array_ptr_h__
#define __retin_auto_array_ptr_h__

#include "core.h"

#include <vector>

namespace retin {

template<class T>
class	auto_array_ptr
{
	T*	ptr;
	size_t	nSize;
	bool	bOwnData;
public:
	auto_array_ptr(size_t nSize=0) : nSize(nSize)
	{
		if (nSize > 0)
		{
                    try {
                        ptr = new T [nSize];
                        if (!ptr)
                            retinThrowException1("Bad alloc size %lu",(unsigned long)nSize);
                    }
                    catch (const std::bad_alloc& ex) {
                        retinThrowException1("Bad alloc size %lu",(unsigned long)nSize);
                    }
                    memset (ptr,0,nSize*sizeof(T));
                    bOwnData = true;
		}
		else
		{
                    ptr = NULL;
                    bOwnData = false;
		}
	}
	auto_array_ptr(T* ptr,size_t nSize,bool bOwnData=false) : ptr(ptr),nSize(nSize),bOwnData(bOwnData)
	{
	}
	auto_array_ptr(const auto_array_ptr<T>& x) : ptr(NULL),nSize(0),bOwnData(false) { copy(x); }
    auto_array_ptr(const std::vector<T>& x) : ptr(NULL),nSize(0),bOwnData(false) { copy(&x[0],x.size()); }
	~auto_array_ptr()
	{
		clear();
	}
	
	T*		get(size_t i=0) { return ptr+i; }
	
	T* const	get(size_t i=0) const { return ptr+i; }
	
	T*		release() { bOwnData = false; return ptr; }
	
	size_t		size() const { return nSize; }
    
    size_t          memoryUsage() const { return nSize*sizeof(T); }
    
	void		clear ()
	{
		if (ptr && bOwnData)
			delete [] ptr;
		nSize = 0;
		ptr = NULL;
                nSize = 0;
		bOwnData = false;
	}
    
    void            take(T* p,size_t n) {
            clear();
            if (p) {
                nSize = n;
                bOwnData = true;
                ptr = p;
            }
        }
        
	void		resize(size_t n)
	{
		 if (nSize == n)
			return;
		clear ();
		nSize = n;
		if (nSize > 0)
		{
			bOwnData = true;
                        try {
                            ptr = new T [nSize];
                            if (!ptr)
                                retinThrowException1("Bad alloc size %lu",(unsigned long)nSize);
                            memset (ptr,0,nSize*sizeof(T));
                        }
                        catch (const std::bad_alloc& ex) {
                            retinThrowException1("Bad alloc size %lu",(unsigned long)nSize);
                        }
		}
	}
	void		copy(const auto_array_ptr<T>& x)
	{
            copy (x.ptr,x.nSize);
	}
        void            copy(const T* p,size_t n) {
            if (!p) {
                clear();
                return;
            }
            if (!ptr || nSize != n) {
                if (ptr && bOwnData)
                        delete [] ptr;
                nSize = n;
                bOwnData = true;
                ptr = new T [nSize];
            }
            memcpy (ptr,p,nSize*sizeof(T));
        }
	void		set (T* p)
	{
		clear();
		ptr = p;
	}
	void		set (T* p,size_t n)
	{
		clear();
		nSize = n;
		ptr = p;
	}
	operator bool() const { return ptr != NULL; }
	void		operator= (const auto_array_ptr<T>& x) { copy(x); }
	T&		operator[] (size_t i) { return ptr[i]; }
	const T&	operator[] (size_t i) const { return ptr[i]; }

};

}

#endif
