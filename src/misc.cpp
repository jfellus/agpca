/*
 * misc.cpp
 *
 *  Created on: 13 nov. 2013
 *      Author: jfellus
 */




///////////
// CONNECTIVITY




void dump_degrees() {
	for(int i=0; i<N; i++) fappend("data/stats/degrees.txt", fmt("%u\n", node[i].nbNeighbors));
}


void create_network() {
	node = new Node[N];
	if(TOPO=="full" || TOPO=="BC") {
		DBG("Fully connected network");
		for(int i=0; i<N; i++) {	node[i].nbNeighbors = N-1; node[i].neighbors = NULL;	}
		NB_EDGES = N*(N-1);
	} else if(TOPO=="BA") {
		DBG("BA-model network");
		DBG("  max degree = " << NB_NEIGHBORS);
		// BA-model (preferential attachment)
		for(int i=0; i<N; i++) { node[i].neighbors = new int[NB_NEIGHBORS]; node[i].nbNeighbors = 0; }
		node[0].connect(1);

		for(int i=2; i<N; i++) {
			for(int j=0; j<i-1; j++) {
				if(node[j].nbNeighbors >= NB_NEIGHBORS) continue;
				if(node[j].nbNeighbors == 0) continue;
				int e = rand()%NB_EDGES;
				if(e<=node[j].nbNeighbors) node[i].connect(j);
			}
			node[i].connect(i-1);
		}

		dump_degrees();
	} else if(TOPO=="MW") {
		node[0].neighbors = new int[N-1]; node[0].nbNeighbors = N-1;
		for(int i=1; i<N; i++) {
			node[i].neighbors = new int[1];
			node[i].nbNeighbors = 1;
			node[i].neighbors[0] = 0;
			node[0].neighbors[i-1] = i;
		}
	} else if(TOPO=="ST") {
		for(int i=0; i<N; i++) {
			int nn = 0;
			if(i!=0) nn ++;
			if(2*i + 1 < N) nn++;
			if(2*i + 2 < N) nn++;
			node[i].neighbors = new int[nn];
			node[i].nbNeighbors = nn;
			if(i!=0) node[i].neighbors[node[i].nbNeighbors-1] = (i-1)/2;
			if(2*i + 1 < N) node[i].neighbors[0] = 2*i + 1;
			if(2*i + 2 < N) node[i].neighbors[1] = 2*i + 2;
		}
	}
}

void generateData(const char* s) {
	X.randGaussian(D,n,p);
}




void init(const char* datafile) {
	DBG("INIT");
	DBGV(NBTHREADS);
	system("rm -rf data/*");
	system("rm -rf plots/*");

	if(datafile[0]=='*') {
		generateData(&datafile[1]);
	}
	else X.load(datafile);
	if(LIMIT_NDATA!=-1 && X.height > LIMIT_NDATA) X.height = LIMIT_NDATA;
	n = X.height;
	D = X.width;

	global_mu.create(D,1);

	FULL_GRAM.create(D,D);
	FULL_COVARIANCE_th.create(D,D);


	create_network();

	int ndo = 0;
	for(int i=0; i<N-1; i++) {
		node[i].init(X, ndo, n/N);
		ndo += n/N;
	}
	node[N-1].init(X,ndo, n-ndo);

//	mat_plotall.create((D+1)*(int)(sqrt(N)+1),(D+1)*(int)(sqrt(N)+1));

	DBGV(LIMIT_NDATA);
	DBGV(N);
	DBGV(D);
	DBGV(n);
	if(USE_ENERGY) {DBGV(KEEP_ENERGY);}
	else {DBGV(q);}
}




void deinit() {
	delete[] node;
}
