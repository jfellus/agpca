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
 * \file work_set.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __retin_work_set_h__
#define __retin_work_set_h__

#include "work_pool.h"
#include "work_runnable.h"

namespace retin {

class   work_set
{
public:
    
        class conqueror {
        public:
            virtual ~conqueror() { }
            virtual void    conquer (size_t i0,size_t n) = 0;
        };

private:
        work_pool*      pool;
        size_t          counter;

	template<typename F>
	class func_task : public work_task {
	    F func;
	public:
	    func_task(const F& func,size_t& counter) : work_task(counter),func(func) { }
	    virtual void	run() { func(); }
	};
        
        class class_task : public work_task {
	    work_runnable& runnable;
	public:
	    class_task(work_runnable& runnable,size_t& counter) : work_task(counter),runnable(runnable) { }
	    virtual void	run() { runnable.run(); }
	};
     
        class conqueror_task : public work_task {
            conqueror& c;
            size_t i0,n;
        public:
            conqueror_task(conqueror& c,size_t i0,size_t n,size_t& counter) : work_task(counter),c(c),i0(i0),n(n) {}
            virtual void	run() { c.conquer(i0,n); }
        };       
        
public:

        #ifdef RETIN_ENABLE_BOOST
        void divide_and_conquer(conqueror& c,size_t i0,size_t size) {
            if (!pool) {
                c.conquer(i0,size);
                return;
            }
            size_t j = 0;
            size_t taskCount = 2 * threadCount();
            size_t m = size / taskCount;
            if ( (size%taskCount) != 0)
                m ++;
            while( j < size ) {
                if ( (j+m) > size )
                    m = size - j;
                if (m > 0) {
                    pool->queue.push(new conqueror_task(c,i0+j,m,counter)); 
                }
                j += m;
            }
        }
        #else
        void divide_and_conquer(conqueror& c,size_t i0,size_t size) {
            c.conquer(i0,size);
        }
        #endif

        work_set(work_pool* pool) : pool(pool),counter(0) {
        }
virtual ~work_set() {
        }

        work_pool* workpool() const {
            return pool;
        }
        size_t  threadCount() const {
            if (!pool) return 1;
            return pool->threadCount();
        }
        size_t  size() const {
            return counter;
        }


        #ifdef RETIN_ENABLE_BOOST
	template<typename F>
	void	push (const F& func) {
            if (pool) {
                pool->queue.push(new func_task<F>(func,counter));
            }
            else
                func_task<F>(func,counter).run();
	}
        
        void    pushRunnable(work_runnable & runnable) {
            if(pool) {
                pool->queue.push(new class_task(runnable,counter));
            }
            else {
                class_task(runnable,counter).run();
            }
        }
        
	void	join() {
            if (!pool || counter == 0)
                return;
            pool->run();
            pool->queue.wait(counter);
	}
        #else
	template<typename F>
	void	push (const F& func) {
            func_task<F>(func,counter).run();
	}
        void	pushRunnable (work_runnable & runnable) {
            class_task(runnable,counter).run();
	}
	void	join() {
            std::cout << "No Multithreading !" << std::endl;
	}
        #endif
};

}

#endif
