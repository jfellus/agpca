/*
 * multithread.cpp
 *
 *  Created on: 22 nov. 2013
 *      Author: jfellus
 */


#include "multithread.h"
#include "utils.h"


int get_nb_cores() {
#ifdef MONOTHREAD
	return 1;
#else
#ifdef _SC_NPROCESSORS_ONLN
  return sysconf(_SC_NPROCESSORS_ONLN);
#else
  	  return 1;
#endif
#endif
}

int NBTHREADS = get_nb_cores();


int MULTITHREAD_NB_DONE = 0;
int MULTITHREAD_NB_PROCESSED = 0;
int MULTITHREAD_TOTAL_TOPROCESS = 0;
pthread_cond_t MULTITHREAD_COND = PTHREAD_COND_INITIALIZER;
pthread_mutex_t MULTITHREAD_MUT = PTHREAD_MUTEX_INITIALIZER;

typedef struct {void(*f)(int); int first,nb;} _pthread_do_multithread_data;
typedef struct {void(*f)(int,int); int first,nb;} _pthread_do_multithread_range_data;

void* _pthread_do_multithread_range(void* d) {
	_pthread_do_multithread_range_data* data = (_pthread_do_multithread_range_data*)d;

	data->f(data->first,data->nb);

//	delete data;

	pthread_mutex_lock(&MULTITHREAD_MUT);
	MULTITHREAD_NB_DONE++;
	pthread_cond_signal(&MULTITHREAD_COND);
	pthread_mutex_unlock(&MULTITHREAD_MUT);

	return 0;
}
void* _pthread_do_multithread(void* d) {
	_pthread_do_multithread_data* data = (_pthread_do_multithread_data*)d;

	for(int i=0; i<data->nb; i++) {
		data->f(i+data->first);
		pthread_mutex_lock(&MULTITHREAD_MUT);
		MULTITHREAD_NB_PROCESSED++;
		DBG_PERCENT((float)MULTITHREAD_NB_PROCESSED/MULTITHREAD_TOTAL_TOPROCESS);
		pthread_mutex_unlock(&MULTITHREAD_MUT);
	}

//	delete data;

	pthread_mutex_lock(&MULTITHREAD_MUT);
	MULTITHREAD_NB_DONE++;
	pthread_cond_signal(&MULTITHREAD_COND);

	pthread_mutex_unlock(&MULTITHREAD_MUT);

	return 0;
}


void do_multithread(void (*f)(int i), int n) {
	DBG("__");
	int nbthreads = NBTHREADS;
	if(nbthreads>n) nbthreads = n;

	pthread_t t[nbthreads];
	int first = 0;
	int nb = n/nbthreads;

	MULTITHREAD_NB_DONE = 0;
	MULTITHREAD_TOTAL_TOPROCESS = n;
	MULTITHREAD_NB_PROCESSED = 0;

	for(int i=0; i<nbthreads; i++) {
		if(first+nb > n) nb = n - first;
		_pthread_do_multithread_data* d = new _pthread_do_multithread_data;
		d->f = f; d->first = first; d->nb = nb;
		pthread_create(&t[i], 0, _pthread_do_multithread, d);
		first += nb;
	}

	pthread_mutex_lock(&MULTITHREAD_MUT);
	while(MULTITHREAD_NB_DONE<nbthreads) {
		pthread_cond_wait(&MULTITHREAD_COND, &MULTITHREAD_MUT);
	}
	pthread_mutex_unlock(&MULTITHREAD_MUT);
}


void do_multithread(void (*f)(int start,int end), int n) {
	int nbthreads = NBTHREADS;
	if(nbthreads>n) nbthreads = n;

	pthread_t t[nbthreads];
	int first = 0;
	int nb = n/nbthreads;

	MULTITHREAD_NB_DONE = 0;

	for(int i=0; i<nbthreads; i++) {
		if(first+nb > n) nb = n - first;
		_pthread_do_multithread_range_data* d = new _pthread_do_multithread_range_data;
		d->f = f; d->first = first; d->nb = nb;
		pthread_create(&t[i], 0, _pthread_do_multithread_range, d);
		usleep(1000);
	}

	pthread_mutex_lock(&MULTITHREAD_MUT);
	while(MULTITHREAD_NB_DONE<nbthreads) {
		pthread_cond_wait(&MULTITHREAD_COND, &MULTITHREAD_MUT);
	}
	pthread_mutex_unlock(&MULTITHREAD_MUT);
}

void CRITICAL_BEGIN() {
	pthread_mutex_lock(&MULTITHREAD_MUT);
}

void CRITICAL_END() {
	pthread_mutex_unlock(&MULTITHREAD_MUT);
}

int get_thread_id() {
	return (int) pthread_self();
}
