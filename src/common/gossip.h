/*
 * gossip.h
 *
 *  Created on: 9 nov. 2013
 *      Author: jfellus
 */

#ifndef GOSSIP_H_
#define GOSSIP_H_

extern int N;

int gossip_choose_sender() {
	return rand()%N;
}

int gossip_choose_receiver(int sender) {
	int r = rand()%(N-1);
	if(r>=sender) r++;
	return r;
}





#endif /* GOSSIP_H_ */
