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
 * \file work_queue.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __retin_work_queue_h__
#define __retin_work_queue_h__

#include "work_task.h"

#include <boost/thread.hpp>
#include <queue>

namespace retin {

class work_queue {

public:

protected:
    std::queue<work_task*>		the_queue;
    mutable boost::mutex	the_mutex;
    boost::condition_variable   work_task_pushed;
    boost::condition_variable   work_task_done;

public:

    void push(work_task* t)
    {
        boost::mutex::scoped_lock lock(the_mutex);
        t->counter ++;
        //std::cout << "push() counter = " << t->counter << std::endl;
        the_queue.push(t);
        lock.unlock();
        work_task_pushed.notify_one();
    }
    size_t size() const
    {
	return the_queue.size();
    }
    bool empty() const
    {
        boost::mutex::scoped_lock lock(the_mutex);
        return the_queue.empty();
    }

    work_task* wait_and_pop()
    {
        boost::mutex::scoped_lock lock(the_mutex);
        while(the_queue.empty())
        {
            work_task_pushed.wait(lock);
        }
        work_task* t = the_queue.front();
        the_queue.pop();
        return t;
    }

    void notify_task_done(work_task* t)
    {
        boost::mutex::scoped_lock lock(the_mutex);
        //std::cout << "pop() counter = " << t->counter << std::endl;
        t->counter --;
        lock.unlock();
        t->done();
        work_task_done.notify_all();
    }

    void wait(size_t& counter)
    {
        boost::mutex::scoped_lock lock(the_mutex);
        while (counter) {
            work_task_done.wait(lock);
            //std::cout << "wait() counter = " << counter << std::endl;
        }
    }

};



}

#endif
