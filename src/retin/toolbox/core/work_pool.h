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
 * \file work_pool.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __retin_work_pool_h__
#define __retin_work_pool_h__


#ifdef RETIN_ENABLE_BOOST

#include "work_queue.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace retin {

    
class work_pool {

friend class work_set;

	class final_task : public work_task {
        public:
            final_task(size_t& counter) : work_task(counter) { }
	    virtual void	run() { }
            virtual bool	stop() { return true; }
	};

	class thread {
	    work_queue&	queue;
	public:
	    thread(work_queue& queue) : queue(queue) {
	    }
	    void	operator()() {
		while (true) {
		    work_task* t = queue.wait_and_pop();
                    if (t->stop()) {
                        delete t;
                        break;
                    }
                    t->run();
                    queue.notify_task_done(t);
                    delete t;
		}
	    }
	};

	size_t	nThreadCount;
	work_queue queue;
	boost::shared_ptr<boost::thread_group> pool;

public:
	work_pool(size_t n=0) : nThreadCount(0) {
            setThreadCount(n);
	}
	~work_pool() {
            destroy();
	}
	size_t  threadCount() {
	    return nThreadCount;
	}
	size_t	taskCount() {
	    return queue.size();
	}
        void    setThreadCount(size_t n) {
            size_t cpuCount = boost::thread::hardware_concurrency();
            if (n == 0 || n > cpuCount)
                n = cpuCount;
            if (n != nThreadCount) {
                destroy();
                nThreadCount = n;
            }
        }
	void	run() {
            if (!pool) {
                pool = boost::make_shared<boost::thread_group>();
                for (size_t i=0;i<nThreadCount;i++) {
                    pool->create_thread(thread(queue));
                }
            }
	}
        void    destroy() {
            if (!pool)
                return;
            size_t counter = 0;
	    for (size_t i=0;i<pool->size();i++)
                queue.push(new final_task(counter));
            pool->join_all();
            pool.reset();
        }
};

}

#else

#include "work_task.h"

namespace retin {

    
class work_pool {

friend class work_set;

public:
	work_pool(size_t n=0) {
	}
	~work_pool() {
            destroy();
	}
	size_t  threadCount() {
	    return 1;
	}
	size_t	taskCount() {
	    return 1;
	}
        void    setThreadCount(size_t n) {
        }
	void	run() {
	}
        void    destroy() {
        }
};

}

#endif

#endif
