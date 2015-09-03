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
 * \file arrays.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __retin_arrays_h__
#define __retin_arrays_h__

#include "core.h"

#include <vector>
#include <algorithm>

namespace retin {

	template<class Array>
	struct index_compare {
		const Array& array;
		index_compare(const Array& array) : array(array) { }
		bool	operator() (const size_t a,const size_t b) { return array[a] < array[b]; }
	};
	template<class T>
	void	isort_inc (std::vector<size_t>& indexes,const std::vector<T>& list) {
		indexes.resize(list.size());
		for (size_t i=0;i<indexes.size();i++) indexes[i] = i;
		sort (indexes.begin(),indexes.end(),index_compare< std::vector<T> >(list));
	}
	template<class Array>
	struct index_rcompare {
		const Array& array;
		index_rcompare(const Array& array) : array(array) { }
		bool	operator() (const size_t a,const size_t b) { return array[a] > array[b]; }
	};
	template<class T>
	void	isort_dec (std::vector<size_t>& indexes,const std::vector<T>& list) {
		indexes.resize(list.size());
		for (size_t i=0;i<indexes.size();i++) indexes[i] = i;
		sort (indexes.begin(),indexes.end(),index_rcompare< std::vector<T> >(list));
	}
	template<class T>
	void	iperm (std::vector<T>& list,const std::vector<size_t>& indexes) {
		std::vector<T> temp(list);
		for (size_t i=0;i<indexes.size();i++) list[i] = temp[indexes[i]];
	}


	template<class K,class T>
	bool	compare_first (const std::pair<K,T>& a,const std::pair<K,T>& b) {
		return a.first < b.first;
	}
	template<class K,class T>
	void	ksort (std::vector< std::pair<K,T> >& array) {
		sort (array.begin(),array.end(),compare_first<K,T>);
	}

	template<class K,class T>
	bool	compare_second_dec (const std::pair<K,T>& a,const std::pair<K,T>& b) {
		return a.second > b.second;
	}
	//! Tri les valeurs par ordre croissant.
	template<class K,class T>
	void	vsort_dec (std::vector< std::pair<K,T> >& array) {
		sort (array.begin(),array.end(),compare_second_dec<K,T>);
	}
	//! Tri les valeurs par ordre croissant.
	template<class K,class T>
	void	partial_vsort_dec (std::vector< std::pair<K,T> >& array,size_t k) {
		partial_sort (array.begin(),array.begin()+k,array.end(),compare_second_dec<K,T>);
	}

	template<class K,class T>
	bool	compare_second_inc (const std::pair<K,T>& a,const std::pair<K,T>& b) {
		return a.second < b.second;
	}
	//! Tri les valeurs par ordre décroissant.
	template<class K,class T>
	void	vsort_inc (std::vector< std::pair<K,T> >& array) {
		sort (array.begin(),array.end(),compare_second_inc<K,T>);
	}
	//! Tri les valeurs par ordre décroissant.
	template<class K,class T>
	void	partial_vsort_inc (std::vector< std::pair<K,T> >& array,size_t k) {
		partial_sort (array.begin(),array.begin()+k,array.end(),compare_second_inc<K,T>);
	}

}

inline std::ostream& operator << (std::ostream& os, std::vector<float>& v) {
	for (size_t i=0;i<v.size();i++)
		os << v[i] << ", ";
	return os;
}


#endif
