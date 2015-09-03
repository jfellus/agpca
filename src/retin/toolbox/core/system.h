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
 * \file system.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __retin_system_h__
#define __retin_system_h__

#include "core.h"

#include <boost/system/config.hpp>

#ifdef BOOST_POSIX_API

#include <sys/stat.h>
#include <stdint.h>
#include <sys/time.h>
#include <sys/resource.h>
namespace retin {
    inline time_t	getFileChangeTime (const char* pcFileName)
    {
	    struct stat info;
	    if (stat(pcFileName,&info) != 0)
		    throw std::runtime_error ("Erreur analyse date fichier "+std::string(pcFileName));
	    return info.st_mtime;
    }

    inline double getcputime(void)
    {
	    struct timeval tim;
	    struct rusage ru;
	    getrusage(RUSAGE_SELF, &ru);
	    tim = ru.ru_utime;
	    double t=(double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;
	    /*tim = ru.ru_stime;
	    t += (double)tim.tv_sec + (double)tim.tv_usec / 1000000.0;*/
	    return t;
    }
    inline double getdaytime()
    {
	    struct timeval tv;
	    gettimeofday(&tv,NULL);
	    return double(tv.tv_sec)+1E-6*double(tv.tv_usec);
    }
}
#else

namespace retin {
	inline time_t	getFileChangeTime (const char* pcFileName)
	{
		retinThrowException("Not implemented");
	}

	inline double getcputime(void)
	{
		retinThrowException("Not implemented");
	}
}
#endif

#endif
