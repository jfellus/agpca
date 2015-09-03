/*
 * multithread.h
 *
 *  Created on: 8 nov. 2013
 *      Author: jfellus
 */

#ifndef MULTITHREAD_H_
#define MULTITHREAD_H_

#include <pthread.h>
#include <time.h>
#include <unistd.h>

int get_nb_cores();

extern int NBTHREADS;


extern int MULTITHREAD_NB_DONE;
extern pthread_cond_t MULTITHREAD_COND;
extern pthread_mutex_t MULTITHREAD_MUT;

void* _pthread_do_multithread_range(void* d);
void* _pthread_do_multithread(void* d);


void do_multithread(void (*f)(int i), int n);


void do_multithread(void (*f)(int start,int end), int n);



#define __multithread__(func)  \
		void _do_##func (int); void func(int N) {do_multithread(_do_##func, N);}	\
		void _do_##func

void CRITICAL_BEGIN();

void CRITICAL_END();

int get_thread_id();

#endif /* MULTITHREAD_H_ */
