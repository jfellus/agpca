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
 * \file cache_queue.h
 * \author Philippe H. Gosselin
 * \version 4.0
 */

#ifndef __retin_cache_queue_h__
#define __retin_cache_queue_h__

namespace retin {

/*!
 * Classe pour tenir à jour une liste donnée (type Data) triée en fonction du nombre d'accès.
 *
 * Toutes les fonctions ont une complexité en O(1)
 *
 * offer(index) :
 *	ajoute une nouvelle donnée identifiée par "index" si elle n'existe pas
 *	sinon elle renvoie la donnée actuelle correspondant à "index"
 *	cette donnée est passée en tête
 *
 * getData(index) :
 *	renvoie la donnée correspondant à "index"
 *
 * getData(index,true) :
 *	renvoie la donnée correspondant à "index" et la passe en tête
 *
 * getHead() :
 *	renvoie la donnée la plus souvent demandée
 *
 * getTail() :
 *	renvoie la donnée la moins souvent demandée
 *
 **/
template<class Data>
class cache_queue
{
protected:
	class Node {
	public:
		Node*		prev;
		Node*		next;
		size_t		index;
		Data		data;
		Node*		nextData;
		Node() : prev(NULL),next(NULL),index(0) { }
		Node(size_t index) : prev(NULL),next(NULL),index(index) { }
	};
	size_t		p,count,depth;
	Node*		head;
	Node*		tail;
	Node**		partitions;

	size_t	hash(size_t index) const {
		return (index*41) % p;
	}
	Node*	find (size_t index) const {
		size_t h = hash(index);
		Node* partition = partitions[h];
		while (partition) {
			if (partition->index == index)
				break;
			partition = partition->nextData;
		}
		return partition;
	}
	Node*	insert (size_t index) {
		size_t h = hash(index);
		Node* partition = partitions[h];
		Node* node = new Node(index);
		node->nextData = partition;
		partitions[h] = node;
		count ++;
		return node;
	}
	Node*	del (size_t index) {
		size_t h = hash(index);
		Node* prev = NULL;
		Node* partition = partitions[h];
		while (partition) {
			if (partition->index == index)
				break;
			prev = partition;
			partition = partition->nextData;
		}
		if (!partition)
			return NULL;
		if (prev)
			prev->nextData = partition->nextData;
		else
			partitions[h] = partition->nextData;
		count --;
		return partition;
	}

public:
	cache_queue(size_t p=1000003) : p(p),count(0),depth(0),head(NULL),tail(NULL),partitions(NULL) {
		partitions = new Node*[p];
		for (size_t i=0;i<p;i++)
			partitions[i] = NULL;
	}
	~cache_queue() {
		clear ();
		if (partitions) {
			delete [] partitions;
			partitions = NULL;
		}
	}

	void	clear() {
		Node* finger = head;
		while (finger) {
			Node* next = finger->next;
			delete finger;
			finger = next;
		}
		head = NULL;
		tail = NULL;
		for (size_t i=0;i<p;i++)
			partitions[i] = NULL;
		depth = 0;
		count = 0;
	}
	size_t	size () const {
/*		Node* finger = head;
		size_t i = 0;
		while (finger)
			finger = finger->next;
		return i;*/
		return count;
	}
	size_t	getDepth() const {
		return depth;
	}

	Data*	getData(size_t index) const {
		Node* current = find(index);
		return current?(&current->data):NULL;
	}
	Data*	getData(size_t index,bool bMoveToHead) {
		Node* current = find(index);
		if (current) {
			if (bMoveToHead) {
			    if (head == current)
				    return &current->data;
			    else if (tail == current)
				    tail = current->prev;
			    current->prev->next = current->next;
			    if (current->next)
				    current->next->prev = current->prev;
			    if (head)
				    head->prev = current;
			    else
				    tail = current;
			    current->prev = NULL;
			    current->next = head;
			    head = current;
			}
			return &current->data;
		}
		return NULL;
	}
	Data*	offer(size_t index) {
		Node* current = find(index);
		if (current) {
			if (head == current)
				return &current->data;
			else if (tail == current)
				tail = current->prev;
			current->prev->next = current->next;
			if (current->next)
				current->next->prev = current->prev;
		}
		else {
			current = insert(index);
		}
		if (head)
			head->prev = current;
		else
			tail = current;
		current->prev = NULL;
		current->next = head;
		head = current;
		return &current->data;
	}
	Data*	moveToTail(size_t index) {
		Node* current = find(index);
		if (current) {
		    if (tail == current)
			    return &current->data;
		    else if (head == current)
			    head = current->next;
		    current->next->prev = current->prev;
		    if (current->prev)
			    current->prev->next = current->next;
		    if (tail)
			    tail->next = current;
		    else
			    head = current;
		    current->prev = tail;
		    current->next = NULL;
		    tail = current;
		    return &current->data;
		}
		return NULL;
	}
	bool	remove(size_t index) {
		Node* current = del(index);
		if (current) {
			if (head == current)
				head = current->next;
			else
				current->prev->next = current->next;
			if (tail == current)
				tail = current->prev;
			else
				current->next->prev = current->prev;
			delete current;
			return true;
		}
		return false;
	}
	Data*	getHead() {
		return head?(&head->data):NULL;
	}
	bool	removeHead() {
		if (!head)
			return false;
		return remove (head->index);
	}
	Data*	getTail() {
		return tail?(&tail->data):NULL;
	}
	bool	removeTail() {
		if (!tail)
			return false;
		return remove (tail->index);
	}

	void	showFromHead () const {
		Node* finger = head;
		int i = 0;
		while (finger && i++<100) {
			std::cout << finger->index << ", ";
			finger = finger->next;
		}
		std::cout << "end" << std::endl;
	}
	void	showFromTail () const {
		Node* finger = tail;
		int i = 0;
		while (finger && i++<100) {
			std::cout << finger->index << ", ";
			finger = finger->prev;
		}
		std::cout << "end" << std::endl;
	}

	class iterator {
	friend class cache_queue;
		Node* node;
	public:
		iterator(Node* node) : node(node) { }
		void	operator=(const iterator& ite) { node = ite.node; }
		bool	operator==(const iterator& ite) { return ite.node == node; }
		bool	operator!=(const iterator& ite) { return ite.node != node; }
		Data&	operator*() { return node->data; }
		Data*	operator->() { return &node->data; }
		Data*	data() { return &node->data; }
		void	operator++(int)	{ if (node) node = node->next; }
	};
	iterator    begin() { return iterator(head); }
	iterator    end() { return iterator(NULL); }
	bool	remove (iterator& ite) {
		if (!ite.node)
			return false;
		Node* node = ite.node;
		ite.node = node->next;
		return remove(node->index);
	}

	class riterator {
	friend class cache_queue;
		Node* node;
	public:
		riterator(Node* node) : node(node) { }
		void	operator=(const riterator& ite) { node = ite.node; }
		bool	operator==(const riterator& ite) { return ite.node == node; }
		bool	operator!=(const riterator& ite) { return ite.node != node; }
		Data&	operator*() { return node->data; }
		Data*	operator->() { return &node->data; }
		Data*	data() { return &node->data; }
		void	operator++(int)	{ if (node) node = node->prev; }
	};
	riterator    rbegin() { return riterator(tail); }
	riterator    rend() { return riterator(NULL); }
	bool	remove (riterator& ite) {
		if (!ite.node)
			return false;
		Node* node = ite.node;
		ite.node = node->prev;
		return remove(node->index);
	}

	class const_riterator {
	friend class cache_queue;
		Node* node;
	public:
		const_riterator(Node* node) : node(node) { }
		void	operator=(const const_riterator& ite) { node = ite.node; }
		bool	operator==(const const_riterator& ite) { return ite.node == node; }
		bool	operator!=(const const_riterator& ite) { return ite.node != node; }
		const Data&	operator*() { return node->data; }
		const Data*	operator->() { return &node->data; }
		const Data*	data() { return &node->data; }
		void	operator++(int)	{ if (node) node = node->prev; }
	};
	const_riterator    rbegin() const { return const_riterator(tail); }
	const_riterator    rend() const { return const_riterator(NULL); }
};

}

#endif
