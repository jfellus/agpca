/*
 * plot.h
 *
 *  Created on: 8 nov. 2013
 *      Author: jfellus
 */

#ifndef PLOT_H_
#define PLOT_H_

#include <string>
#include "matrix.h"
#include "multithread.h"
#include "jpg.h"

using namespace std;

inline void system(const string& s) {if(system(s.c_str())) throw "error";}

void plot_histogram(const Matrix& m, const string& file, float min=0, float max = 100) {
	string ff = fmt("tmp%d.txt",get_thread_id());
	m.save(ff);
	system(fmt("python scripts/plot_histogram.py %s %s %f %f", ff.c_str(), file.c_str(), min, max));
	unlink(ff.c_str());
}

void plot_img(const Matrix& m, const string& file, float min=0, float max=100) {
	FILE* f = fopen(file.c_str(), "w");
	write_jpeg(f,m,m.width,m.height,min,max);
	fclose(f);
}


void plot(const Matrix& m, const string& file, float min = 0, float max = 100) {
	if(m.height==1) plot_histogram(m, file, min, max);
	else plot_img(m, file, min, max);
}

void plot3d(const Matrix& m, const string& file, float min = 0, float max = 100) {
	string ff = fmt("tmp%d.txt",get_thread_id());
	m.save(ff);
	system(fmt("python scripts/plot_histogram3d.py %s %s %f %f", ff.c_str(), file.c_str(), min, max));
	unlink(ff.c_str());
}

#endif /* PLOT_H_ */
