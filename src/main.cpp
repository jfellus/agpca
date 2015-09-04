/*
 * main.cpp
 *
 *  Created on: 8 nov. 2013
 *      Author: jfellus
 */
#include "common/config.h"
#define __$(x) DBG(#x); x.dump()
#define DDD(x) DBG(#x); x;
////////////
// PARAMS //
////////////

string DATASET = get_config_string("DATASET", "/home/jfellus/Documents/These/prog/bases/mnist/train.idx3-ubyte");
int N = get_config("N", 100);
int T_MAX = get_config("TMAX", N*40);

bool LATE_PCA = get_config("LATE_PCA", 0);

bool USE_ENERGY = get_config("USE_ENERGY", 0);
float KEEP_ENERGY = get_config("KEEP_ENERGY", 0.99);
int q = get_config("NB_PROJS", 10);


int p = get_config("P", 50);

int LIMIT_NDATA = get_config("LIMIT_NDATA", -1);
string TOPO = get_config_string("TOPO", "full");
int NB_NEIGHBORS = get_config("NB_NEIGHBORS", -1);
int n = get_config("n", 10000);
int D = get_config("D", 200);


bool DEBUG = get_config("DEBUG", true);



//////////
// DATA //
//////////

Matrix X;

Matrix mu;
Matrix G;
Matrix FULL_COVARIANCE_th;

int t = 0;
int last_sender=-1;
int last_receiver=-1;

int last_nb_projs_sent;
int	last_msg_size;
long total_msg_size = 0;

//Matrix mat_plotall;
int NB_EDGES = 0;


//////////////////////

FILE* f_E;


void compute_errors();
void flushall() {fflush(f_E);}


#include "algebra.cpp"
int __node_last_id = 0;



/////////////////////////////



// OK MATLAB !!!!
void QR_algo_sum(Matrix& U1, Matrix &L1, Matrix& U2, Matrix& L2) {
	Matrix Q(U1.width, U1.height);
	Matrix r(L1.width,1);
	Matrix M(Q.width,Q.height);
	Q = U1;
	for(int i=0; i<10; i++) {
		DBGV(i);
		M.clear();
		M.ULUtQ(U1,L1,Q);
		M.ULUtQ(U2,L2,Q);
		M.qr(Q,r);
	}
	U1 = Q;
	L1 = r;
}


/** Eigendecompose the compound matrix A=(U1L1U1')+(U2L2U2') into [U,L] with rank(A) eigenvectors */
void eig_of_sum(Matrix& U, Matrix& L, Matrix& U1, Matrix& L1, Matrix& U2, Matrix& L2) {
	Matrix A = U1; A.diag_mul_sqrt(L1);
	Matrix B = U2; B.diag_mul_sqrt(L2);

	Matrix X; X.make_h_block(A,B);
	Matrix Xt = X.transpose();
	DBG(Xt.height << "x" << Xt.width);
	Xt.PCA(U,L,L1.width+L2.width);
}

/** Eigendecompose the compound matrix A=(U1L1U1')+(U2L2U2') into [U1,L1] with rank(U1) eigenvectors */
void eig_of_sum(Matrix& U1, Matrix& L1, Matrix& U2, Matrix& L2) {
	Matrix A = U1; A.diag_mul_sqrt(L1);
	Matrix B = U2; B.diag_mul_sqrt(L2);

	Matrix X; X.make_h_block(A,B);
	Matrix Xt = X.transpose();
	Xt.PCA(U1,L1,L1.width);
}



//////////
// NODE //
//////////

class Node;
Node* node;

class Node {
public:
	int id;

	Matrix X;
	int n;

	Matrix outU,outL;
	Matrix mu;
	Matrix C;
	Matrix C_rec;
	Matrix U,L;

	float nb_projs_kept;

	float w;

	Node() {this->id = __node_last_id++;n=0;nb_projs_kept=0;w=0;neighbors=NULL;nbNeighbors=0;}
	void init(Matrix& X, int first, int n) {
		this->n = n;
		this->X.create_ref(X.get_row(first),D,n);
	}

	~Node() {}


	// 1) LOCAL PCA

	void local_pca() {
		X.mean_row(mu);
		X.PCA(U,L,q);
		outU = U;outL = L;
	}


	// 2) GOSSIP

	void init_gossip() {
		w = n; 		mu *= w;	L *= w;
	}

	void send(Node& node) {
		w /= 2; 	mu /= 2;	L /= 2;
		node.receive(U,L,mu,w);
	}

	void receive(Matrix& U2, Matrix& L2, Matrix& mu2, float w2) {
		eig_of_sum(U,L,U2,L2);
		w += w2;	mu += mu2;
	}

	void compute_estimate() {
		outL = L; outL /= w;
		Matrix mu2(1,D); for(int i=0; i<D; i++) mu2[i] = mu[i]/w;
		Matrix one(1,1); one.data[0] = 1;
		outU = U;
		eig_of_sum(outU, outL, mu2, one);
	}


	// 3) Project

	void project(Matrix& X, Matrix &compact) {
		compact.project(X, outU);
	}

	// CONNECTIVITY
	int* neighbors;
	int nbNeighbors;

	bool is_network_complete() {return neighbors==NULL;}
	int get_neighbors_count() {return nbNeighbors;}

	void connect(int neighbor) {
		node[neighbor].neighbors[node[neighbor].nbNeighbors++] = id;
		neighbors[nbNeighbors++] = neighbor;
		NB_EDGES += 2;
	}

	int gossip_choose_receiver() {
		if(is_network_complete()) {
			int r = rand()%(N-1);
			if(r>=id) r++;
			return r;
		} else {
			int r = rand()%get_neighbors_count();
			return neighbors[r];
		}
	}

};


#include "misc.cpp"

/////////////////////

__multithread__(local_pca) (int i) {	node[i].local_pca(); }



void iteration() {
	last_sender = gossip_choose_sender();
	last_receiver = node[last_sender].gossip_choose_receiver();

	node[last_sender].send(node[last_receiver]);
}


//////////
// DUMP //
//////////


void plot_Lmean() {
	Matrix Lmean(D,1);
	for(int i=0; i<N; i++) Lmean += node[i].L;
	Lmean.normalize_l1();
	plot(Lmean, fmt("plots/spectrum_%08d.jpg",t), 0,1);
}


void compute_errors() {
	if(last_sender!=-1) node[last_sender].compute_estimate();
	if(last_receiver!=-1) node[last_receiver].compute_estimate();

	if(last_receiver>=0 && last_receiver<1) {
		if(DEBUG) {
			Matrix compact;
			Matrix Xrec;
			Matrix Crec;

			compact.project(X, node[last_receiver].outU);
			Xrec.unproject(compact, node[last_receiver].outU);
			Crec.covariance(Xrec);

			float norm = G.n2();
			float ERR = G.l2(Crec)/norm;
			DBGV(ERR);

			fappend(f_E, fmt("%u %f\n", t, ERR*100));
			fflush(f_E);
		}
		else {
			node[0].outU.save(fmt("data/U_%08d.fvec", t));
			node[0].outL.save(fmt("data/L_%08d.fvec", t));
		}
	}
	//misc_errors();
}

Matrix Y;
Matrix K;
Matrix Xrec;

void readmat(Matrix& m, const char* s) {
	FILE* f = fopen(s, "r");
	if(!f) DBG("ERROR " << s);
	for(int i=0; i<m.width*m.height; i++) fscanf(f, "%f", &m.data[i]);
	fclose(f);
}

void a() {
	Matrix X(20,9);
	Matrix Y(20,9);
	readmat(X, "X.txt");
	readmat(Y, "Y.txt");

DDD(	Matrix CX; CX.correlation(X);)
DDD(	Matrix CY; CY.correlation(Y);)
DDD(	Matrix C = CX+CY;)



	Matrix U1,L1;
	Matrix U2,L2;
DDD(	X.PCA(U1,L1,9); )
DDD(	Y.PCA(U2,L2,9); )

	Matrix rec1; rec1.reconstruct(U1,L1);
	Matrix rec2; rec2.reconstruct(U2,L2);


	Matrix U,L;
	eig_of_sum(U,L,U1,L1,U2,L2);

	__$(U);
	__$(L);

	Matrix rec; rec.reconstruct(U,L);

	__$(C);
	__$(rec);


}








int main(int argc, char **argv) {
	try {


	//	a();return 0;

		DBG_START("Init ");

		if(argc<=1) {
			DBG("Read data from \"" << DATASET << "\"");
			init(DATASET.c_str());
		} else {
			DBG("with : " << argv[1]);
			init(argv[1]);
		}

		sys("rm -rf  data");
		sys("mkdir -p data");

		f_E = fopen("data/E.txt", "w");

		DBG_END();


		DBG(X.height << "x" << X.width);
		DBG_START("Substract mean");
		X.mean_row(mu);
		X -= mu;
		DBG_END();

		if(DEBUG) {
			DBG_START("Compute global PCA");
			Matrix U,L;
			DDD(		X.PCA(U,L, q);				)
			DDD(		G.covariance(X);		)
			DDD(		Y.project(X,U);			)
			DDD( 		Xrec.unproject(Y,U);    )
			DDD(		K.covariance(Xrec);		)
			DBG(G.height << "x" << G.width);
			DBG(K.height << "x" << K.width);
			float norm = G.n2();
			float ERR_pca = G.l2(K)/norm*100;
			DBG("ERROR PCA : " << ERR_pca);
			fappend("data/ERR_pca.txt", fmt("0 %f\n", ERR_pca));
			fappend("data/ERR_pca.txt", fmt("4000 %f\n", ERR_pca));
			DBG_END();
		}

		////////////////////////

		DBG_START("Local PCA");
		for(int i=0; i<N; i++) {DBG_PERCENT(i/N);node[i].local_pca();}
		if(DEBUG) compute_errors();
		DBG_END();

		for(int i=0; i<N; i++) node[i].init_gossip();

		for(t=0; t<T_MAX; t++) {
			DBGV(t);
			iteration();
			compute_errors();
			if(t%50==0) flushall();
		}

		//////////////////////


		//	Ufull.create(D,D); Lfull.create(D,1);
		//	FULL_COVARIANCE.eigendecompose(Ufull, Lfull);
		//	dump_latepca_errors(D);
		fclose(f_E);
		DBG("finished");


	}catch(const char* c) {DBG("ERROR : " << c );}
}
