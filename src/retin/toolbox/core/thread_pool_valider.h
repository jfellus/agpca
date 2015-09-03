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

#ifndef __thread_pool_valider_h__
#define __thread_pool_valider_h__

#include "thread_pool.h"

#include <sys/time.h>
#include <iostream>

using namespace std;

inline void function_work() {
    double sum = 0;
    size_t n = 100000;
    for (size_t k=0;k<1000;k++)
    for (size_t i=0;i<n;i++)
        sum += exp(-double(n-1-i)/double(n));
    cout << "Result: " << sum << endl;
}

class functor_work {
public:
    void operator()() { function_work(); }
};

class runnable_work : public thread_pool::runnable {
public:
    void run() { function_work(); }
};

class conqueror : public thread_pool::conqueror {
    size_t sum;
public:
    conqueror () : sum(0) { }
    void    conquer (size_t i0,size_t n) {
        size_t s;
        for (size_t t=0;t<100000;t++) {
            s = 0;
            for (size_t j=0;j<n;j++) {
                size_t i = i0+j;
                s += i;
            }
        }
        lock();
        sum += s;
        unlock();
    }
    size_t getResult() { return sum; }
};

static inline double getdaytime()
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return double(tv.tv_sec)+1E-6*double(tv.tv_usec);
}

inline void thread_pool_valider(size_t threadCount=0) {

    thread_pool pool(threadCount);
    thread_pool* ppool = &pool; // Put NULL to test single threaded

    cout << "Thread count: " << pool.getThreadCount() << endl;

    cout << "Run pool..." << endl;
    pool.run();
    usleep(100000);

    thread_pool::set workset(ppool);
    cout << "Push function work..." << endl;
    workset.push(function_work);

    cout << "Push functor work..." << endl;
    functor_work functor;
    workset.push(functor);

    cout << "Push runnable work..." << endl;
    runnable_work runnable;
    workset.push_runnable(runnable);

    cout << "Wait for work set to end..." << endl;
    workset.join();
    cout << "Work set finished." << endl;


    cout << "Try with a conqueror..." << endl;
    double start = getdaytime();
    conqueror c;
    size_t n = 1000000;
    c.divide_and_conquer(ppool,0,n);
    cout << "Result = " << c.getResult() << " == " << ((n-1)*n)/2 << " time=" << (getdaytime()-start) << endl;
}

#endif


