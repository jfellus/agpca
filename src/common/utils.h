/*
 * utils.h
 *
 *  Created on: 8 nov. 2013
 *      Author: jfellus
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <string>
#include <stdarg.h>
#include "multithread.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <string>
#include <libgen.h>
#define MIN(a,b) ((a)>(b) ? (a) : (b))

using namespace std;

#define DBG_START(x) \
	do {CRITICAL_BEGIN(); \
	std::cout << x << " ... "; fflush(stdout); \
	CRITICAL_END();} while(0);

#define DBG_END() \
	do {CRITICAL_BEGIN(); \
	std::cout << "ok\n"; fflush(stdout); \
	CRITICAL_END();} while(0);


#define DBG(x) \
	do {CRITICAL_BEGIN(); \
	std::cout << x << "\n"; fflush(stdout); \
	CRITICAL_END();}while(0);

#define DBGV(x) \
	CRITICAL_BEGIN(); \
	std::cout << "" #x " : " << x << "\n"; fflush(stdout); \
	CRITICAL_END();

#define DBGN(x) \
		do { \
			std::cout << x << "\n"; fflush(stdout); \
			std::cout  << "\033[1F\033[2K"; \
		}while(0);

#define DBG_PERCENT(x) DBGN((int)((float)x*1000)/10.0f << "%")

string stringprintf(const char* fmt, ...);
#define fmt(x...) stringprintf(x)

void fappend(const string& filename, const string& line);
void fappend(FILE *f, const string& line);

#define FATAL(x) do { DBG(x); exit(1); } while(0);

#endif /* UTILS_H_ */
