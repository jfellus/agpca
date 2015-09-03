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
#ifndef __retin_string_h__
#define __retin_string_h__

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <sstream>

namespace retin {

static inline std::string string_format(const std::string &fmt, ...) {
	int size = 100;
	std::string str;
	va_list ap;
	while (1) {
		str.resize(size);
		va_start(ap, fmt);
		int n = vsnprintf((char *) str.c_str(), size, fmt.c_str(), ap);
		va_end(ap);
		if (n > -1 && n < size) {
			str.resize(n);
			return str;
		}
		if (n > -1)
			size = n + 1;
		else
			size *= 2;
	}
}

static inline bool string_has_file_ext(const std::string& file_name,
		const std::string& ext) {
	if (file_name.length() < ext.length()) {
		return false;
	}
	return file_name.rfind(ext) == (file_name.length() - ext.length());
}

static inline std::string& string_set_file_ext(std::string& file_name,
		const std::string& ext) {
	size_t idx = file_name.rfind(".");
	if (idx != std::string::npos) {
		file_name.resize(idx);
	}
	return file_name.append(ext);
}

static inline std::vector<std::string> &string_split(const std::string &s,
		char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

static inline std::vector<std::string> string_split(const std::string &s,
		char delim) {
	std::vector < std::string > elems;
	string_split(s, delim, elems);
	return elems;
}

}

#endif
