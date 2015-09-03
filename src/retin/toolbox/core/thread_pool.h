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

#ifndef __thread_pool_h__
#define __thread_pool_h__

#include <unistd.h>
#include <stdexcept>
#include <queue>

#ifdef THREAD_POOL_USE_POSIX_THREADS    
#include <pthread.h>
#endif

#ifndef THREAD_POOL_VERBOSE
#define THREAD_POOL_VERBOSE 0
#endif

#if THREAD_POOL_VERBOSE >= 1
#include <iostream>
#endif

/*! Multithreading computation

Usage:

 * First create a thread_pool & run it:
        thread_pool pool;
        pool.run();
    It is recommended to create a single pool & run it at the begining of the main, and use it in all your program.
    Notes:
     - If you destroy the pool, it will terminate all threads. Otherwise, the threads are always running, and waiting for tasks.
     - When threads are waiting for tasks, they don't use any cpu (as if they were sleeping).
     - Launching and destroying the pool is slow (some millisecond). It is recommended to use work sets (see below).

 * Push some work straight to the pool:
        pool.push( function );  //  function = any function or class with operator() and no arguments (ex: void f(); )
        pool.push_runnable( runnable );  // runnable = any instance of thread_pool::runnable (implements a method void run(); }
        pool.run();
        pool.destroy(); // will wait for the end of threads
    This way, you have to manage the end of tasks, or just destroy the pool (=> better use work sets).

 * Work sets:
        thread_pool::set workset(&pool);
        workset.push( function );  //  function = any function or class with operator() and no arguments (ex: void f(); )
        workset.push_runnable( runnable );  // runnable = any instance of thread_pool::runnable (implements a method void run(); }
        workset.join(); // wait for the end of pushed tasks
    Notes:
     - If your pool is already running, tasks start immediately
     - You can push tasks in several work sets, each one will only wait for its tasks. Note that scheduling order is unpredictable.
     - Tasks in work sets can create work sets etc.
     - Creating a workset without a pool (ie: thread_pool::set workset(NULL);) will acts as if they is no threads at all. 
       This is interesting if you want to check your code in single threaded environment.

 * Conquerors:
    Conquerors automatically split a list of computations.
    First implements conqueror:
        class my_conqueror : public thread_pool::conqueror {
        public:
            void    process (size_t i) {
                // Compute here all operations for index i
            }
        };
    Or, if you need to perform similar operation per batch:
        class my_conqueror : public thread_pool::conqueror {
        public:
            void    conquer (size_t i0,size_t n) {
                // Compute here all operations between indices [i0,i0+n[
            }
        };
    Then, just run (for instance):
        my_conqueror c;
        c.divide_and_conquer(&pool,0,100000);
    If you want to manage the work set yourself, you can:
        set workset(&pool);
        workset.push_conquerors(c,0,100000);
        workset.join();
    Notes:
     - It creates twice more tasks than split (meaning your tasks can have different running times, but not too much)
     - The base class includes locking functions, if you need to merge results, for instance:
        void    conquer (size_t i0,size_t n) {
            // Compute here all operations between indices [i0,i0+n[
            lock();
            // Merge results (this part won't be multi-threaded)
            unlock();
        }
*/
class thread_pool {

public:

    #ifdef THREAD_POOL_USE_POSIX_THREADS
    class condition;
    // Mutex class
    class mutex {
        friend class condition;
        bool locked;
        pthread_mutex_t m;
    public:
        mutex() : locked(false) {
            if (pthread_mutex_init(&m, NULL))
                throw exception("pthread_mutex_init() error");
        }
        ~mutex() {
            if (locked) 
                pthread_mutex_unlock(&m);
            if (pthread_mutex_destroy(&m))
                throw exception("pthread_mutex_destroy() error");
        }
        void lock() {
            if (pthread_mutex_lock(&m))
                throw exception("pthread_mutex_lock() error");
            locked = true;
        }
        void unlock() {
            if (pthread_mutex_unlock(&m))
                throw exception("pthread_mutex_unlock() error");
            locked = false;
        }
    };

    // Condition class
    class condition {
        pthread_cond_t cond;
    public:
        condition() {
            if (pthread_cond_init(&cond,NULL))
                throw exception("pthread_cond_init() error");
        }    
        ~condition() {
            if (pthread_cond_destroy(&cond))
                throw exception("pthread_cond_destroy() error");
        }
        void    notify_one() {
            if (pthread_cond_signal(&cond))
                throw exception("pthread_cond_signal() error");
        }
        void    notify_all() {
            if (pthread_cond_broadcast(&cond))
                throw exception("pthread_cond_broadcast() error");
        }        
        void    wait(mutex& m) {
            if (pthread_cond_wait(&cond,&m.m))
                throw exception("pthread_cond_wait() error");
        }
    };
    #endif

    // Runnable abstract class
    class   runnable {
    public:
        virtual ~runnable() { }
        virtual void    run() = 0;
    };

    // Conqueror abstract class
    class conqueror {
    protected:
    #ifdef THREAD_POOL_USE_POSIX_THREADS    
        mutable mutex locker;
        void    lock() { locker.lock(); }
        void    unlock() { locker.unlock(); }
    #else
        void    lock() { }
        void    unlock() { }
    #endif
    public:
        virtual ~conqueror() { }
        virtual void    process(size_t i) {
            throw exception("conqueror::process() not implemented");
        }
        virtual void    conquer (size_t i0,size_t n) {  
            for (size_t i=0;i<n;i++) {
                process(i0+i);
            }
        }
        void divide_and_conquer (thread_pool* pool,size_t i0,size_t n);
    };

    // Exception class
    class   exception : public std::exception {
        const char* msg;
    public:
        exception(const char* msg) : msg(msg) {}
        virtual const char* what() const throw() {
            return msg; 
        }
    };
protected:


    #ifdef THREAD_POOL_USE_POSIX_THREADS    
    // Task classes
    class   task {
    public:
        virtual ~task() { }
        virtual bool	stop() { return false; }
        virtual void    enqueued() { }
        virtual void	run() = 0;
        virtual void    done() { }
    };
    class final_task : public task {
    public:
        virtual bool	stop() { return true; }
        virtual void	run() { }
    };
    template<class F>
    class function_task : public task {
        F& func;
    public:
        function_task(F& func) : func(func) { }
        virtual void	run() { func(); }
    };
    class runnable_task : public task {
        runnable& func;
    public:
        runnable_task(runnable& func) : func(func) { }
        virtual void	run() { func.run(); }
    };

    // Queue class
    class concurrent_queue {
    protected:
        std::queue<task*>		queue;
        mutable mutex	        m;
        condition               task_pushed;
        condition               task_done;

    public:
        void lock() const {
            m.lock();
        }
        void unlock() const {
            m.unlock();
        }
        void push(task* t) {
            lock();
            queue.push(t);
            t->enqueued();
            task_pushed.notify_one();
            unlock();
        }
        size_t size() const {
            lock();
    	    size_t s = queue.size();
            unlock();
            return s;
        }
        bool empty() const {
            lock();
            bool b = queue.empty();
            unlock();
            return b;
        }
        task* wait_and_pop() {
            lock();
            while(queue.empty()) {
                task_pushed.wait(m);
            }
            task* t = queue.front();
            queue.pop();
            unlock();
            return t;
        }
        void notify_task_done(task* t) {
            lock();
            t->done();
            task_done.notify_all();
            unlock();
        }
        void wait(size_t& counter) {
            lock();
            while (counter) {
                task_done.wait(m);
            }
            unlock();
        }

    };
    // Worker class
    class worker {
        concurrent_queue&  queue;
        size_t      id;
        pthread_t   thread;
        static void* thread_init(void* p) { ((worker*)p)->run(); return NULL; }
    public:
        worker(concurrent_queue& queue,size_t id) : queue(queue),id(id) {
            if (pthread_create(&thread, NULL, thread_init, (void*)this))
                throw exception("pthread_create() error");
        }
        void   run() {
            #if THREAD_POOL_VERBOSE >= 1
            queue.lock();
            std::cout << "Start thread " << id << std::endl;
            queue.unlock();
            #endif
		    while (true) {
		        task* t = queue.wait_and_pop();
                if (t->stop()) {
                    delete t;
                    break;
                }
                t->run();
                queue.notify_task_done(t);
                delete t;
		    }
            #if THREAD_POOL_VERBOSE >= 1
            queue.lock();
            std::cout << "Stop thread " << id << std::endl;
            queue.unlock();
            #endif
            pthread_exit(NULL);
        }
        void join() {
            if (pthread_join(thread, NULL))
                throw exception("pthread_join() error");
        }
    };
    #endif

public:
    // Set class
    class set {
        thread_pool*  pool;
        size_t      counter;
    public:
        set(thread_pool* pool=NULL) : pool(pool),counter(0) { }
        thread_pool* getWorkpool() const {
            return pool;
        }
        size_t  getThreadCount() const {
            if (!pool) return 1;
            return pool->getThreadCount();
        }
        size_t  size() const {
            return counter;
        }
        template<class F>
        void    push (F& func);
        void    push_runnable(runnable& func);
        void    push_runnable_ptr(runnable* func) { push_runnable(*func); }
        void    push_conquerors(conqueror& c,size_t i0,size_t size);
        void    push_conquerors_ptr(conqueror* c,size_t i0,size_t size);
        void    join();
        void    task_enqueued() {
            counter ++;
            #if THREAD_POOL_VERBOSE >= 1
            std::cout << "Task enqueued; counter = " << counter << std::endl;
            #endif
        }
        void    task_done() {
            counter --;
            #if THREAD_POOL_VERBOSE >= 1
            std::cout << "Task done; counter = " << counter << std::endl;
            #endif
        }
    };

protected:

    #ifdef THREAD_POOL_USE_POSIX_THREADS    
    class set_task : public task {
        set& workset;
    public:
        set_task(set& workset) : workset(workset) { }
        virtual void    enqueued() { workset.task_enqueued(); }
        virtual void	run() = 0;
        virtual void    done() { workset.task_done(); } 
    };
    template<class F>
    class function_set_task : public set_task {
        F& func;
    public:
        function_set_task(set& workset,F& func) : set_task(workset),func(func) { }
        virtual void	run() { func(); }
    };
    class runnable_set_task : public set_task {
        runnable& func;
    public:
        runnable_set_task(set& workset,runnable& func) : set_task(workset),func(func) { }
        virtual void	run() { func.run(); }
    };
    class conqueror_set_task : public set_task {
        conqueror& c;
        size_t i0,n;
    public:
        conqueror_set_task(set& workset,conqueror& c,size_t i0,size_t n) : set_task(workset),c(c),i0(i0),n(n) { }
        virtual void	run() { c.conquer(i0,n); }
    };

	concurrent_queue   queue;
    std::vector<worker*> workers;
    #endif

	size_t	threadCount;
    size_t logicalCoreCount;
public:
	thread_pool(size_t n=0) : threadCount(0) {
        #ifdef THREAD_POOL_USE_POSIX_THREADS    
        logicalCoreCount = sysconf(_SC_NPROCESSORS_CONF);
        #else
        logicalCoreCount = 1;
        #endif
        setThreadCount(n);
	}
	~thread_pool() {
        destroy();
	}
	size_t  getThreadCount() const {
	    return threadCount;
	}
    size_t  getLogicalCoreCount() const {
        return logicalCoreCount;
    }
    void    setThreadCount(size_t n) {
        if (n == 0 || n > logicalCoreCount)
            n = logicalCoreCount;
        if (n != threadCount) {
            destroy();
            threadCount = n;
        }
    }
	void	run() {
        #ifdef THREAD_POOL_USE_POSIX_THREADS    
        if (workers.size() == 0) {
            workers.resize(threadCount);
            for (size_t t=0;t<workers.size();t++) {
                workers[t] = new worker(queue,t);
            }
        }
        #endif
	}
    void    destroy() {
        #ifdef THREAD_POOL_USE_POSIX_THREADS    
        if (workers.size() == 0)
            return;
        #if THREAD_POOL_VERBOSE >= 1
        std::cout << "Wait for threads to end..." << std::endl;
        #endif
        for (size_t t=0;t<workers.size();t++) {
            queue.push(new final_task());
        }
        for (size_t t=0;t<workers.size();t++) {
            workers[t]->join();
        }
        for (size_t t=0;t<workers.size();t++) {
            delete workers[t];
        }
        workers.clear();
        #endif
    }
    #ifdef THREAD_POOL_USE_POSIX_THREADS    
    template<class F>
    void    push (F& func) {
        queue.push(new function_task<F>(func));
    }
    template<class F>
    void    push_set (set& workset,F& func) {
        queue.push(new function_set_task<F>(workset,func));
    }
    void    push_runnable(runnable& func) {
        queue.push(new runnable_task(func));
    }    
    void    push_set_runnable(set& workset,runnable& func) {
        queue.push(new runnable_set_task(workset,func));
    }
    void    push_set_conqueror(set& workset,conqueror& c,size_t i0,size_t n) {
        queue.push(new conqueror_set_task(workset,c,i0,n));
    }
    #else
    template<class F>
    void    push (F& func) {
        func();
    }
    template<class F>
    void    push_set (set& workset,F& func) {
        func();
    }
    void    push_runnable(runnable& func) {
        func.run();
    }    
    void    push_set_runnable(set& workset,runnable& func) {
        func.run();
    }
    void    push_set_conqueror(set& workset,conqueror& c,size_t i0,size_t n) {
        c.conquer(i0,n);
    }
    #endif

    void    divide_and_conquer(conqueror& c,size_t i0,size_t size) {
        set workset(this);
        workset.push_conquerors(c,i0,size);
        workset.join();
    }

};

inline void    thread_pool::conqueror::divide_and_conquer(thread_pool* pool,size_t i0,size_t size) {
    set workset(pool);
    workset.push_conquerors(*this,i0,size);
    workset.join();
}

template<class F>
void    thread_pool::set::push (F& func) {
    #ifdef THREAD_POOL_USE_POSIX_THREADS    
    if (pool)
        pool->push_set(*this,func);
    else
    #endif
        func();
}

inline void    thread_pool::set::push_runnable(runnable& func) {
    #ifdef THREAD_POOL_USE_POSIX_THREADS    
    if (pool)
        pool->push_set_runnable(*this,func);
    else
    #endif
        func.run();
}

inline void     thread_pool::set::push_conquerors(conqueror& c,size_t i0,size_t size) {
    if (!pool) {
        c.conquer(i0,size);
        return;
    }
    size_t j = 0;
    size_t taskCount = 2 * getThreadCount();
    size_t m = size / taskCount;
    if ( (size%taskCount) != 0)
        m ++;
    while( j < size ) {
        if ( (j+m) > size )
            m = size - j;
        if (m > 0) {
            pool->push_set_conqueror(*this,c,i0+j,m); 
        }
        j += m;
    }
}

inline void     thread_pool::set::push_conquerors_ptr(conqueror* c,size_t i0,size_t size) {
    push_conquerors(*c,i0,size);
}

inline void     thread_pool::set::join() {
    #ifdef THREAD_POOL_USE_POSIX_THREADS    
    if (!pool || counter == 0)
        return;
    pool->run();
    pool->queue.wait(counter);
    #endif
}

typedef thread_pool::set thread_pool_set;
typedef thread_pool::runnable thread_pool_runnable;
typedef thread_pool::conqueror thread_pool_conqueror;

#endif


