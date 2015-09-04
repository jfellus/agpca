/*
 * config.h
 *
 *  Created on: 3 sept. 2015
 *      Author: jfellus
 */

#ifndef cCONFIG_H_
#define cCONFIG_H_

#define MONOTHREAD

#include "math.h"
#include "multithread.h"
#include "utils.h"
#include "plot.h"
#include "gossip.h"
#include <vector>
#include <string>
#include <libgen.h>

using namespace std;


template <class T> T get_config(const char* what, T default_val) {
	FILE* f = fopen("config.properties", "r");
	char line[512];
	char* c = 0;
	T v = default_val;
	while ( fgets (line , 512 , f) != NULL ) {
	      if(strlen(what)<strlen(line)
	    &&	!strncmp(what, line, strlen(what))
	    && (line[strlen(what)]==' ' || line[strlen(what)]=='=')
	      ) {
	    	  c = line + strlen(what);
	    	  while(*c==' ' || *c=='=') c++;
	    	  v = (T) atof(c);
	      }
	}
	fclose(f);
	return v;
}

string get_config_string(const char* what, string default_val) {
	FILE* f = fopen("config.properties", "r");
	char line[512];
	char* c = 0;
	string s = default_val;
	while ( fgets (line , 512 , f) != NULL ) {
	      if(!strncmp(what, line, strlen(what))) {
	    	  c = line + strlen(what);
	    	  while(*c==' ' || *c=='=') c++;
	    	  if(c[strlen(c)-1]=='\n') c[strlen(c)-1] = 0;
	    	  s = c;
	      }
	}
	fclose(f);
	return s;
}




#endif /* CONFIG_H_ */
