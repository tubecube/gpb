#ifndef GPB_H
#define GPB_H

#include <vector>
#include <armadillo>
#include <random>
#include <fstream>
#include "graph.hpp"
#include "utils.hpp"

using namespace std;
using namespace arma;

class GPB {

public:
    GPB(Graph& g, int k, double Fpshp=0.3, double Fprte=1.0, double Bpshp=0.3, double Bprte=1.0, const char* dir=".", bool sample_deep=true, double alpha=1.0, double beta=0.0):
    	graph(g), K(k), Fpshp(Fpshp), Fprte(Fprte), Bpshp(Bpshp), Bprte(Bprte), save_dir(dir), sample_deep(sample_deep), alpha(alpha), beta(beta) 
	{ 
		Logger::setup_log_dir(save_dir);
		Logger::setup_logfd(save_dir+"/log.txt");
		INFO("GPB: %s graph with %d nodes %d edges\n", (graph.is_directed() ? "directed" : "undirected"), graph.n_nodes(), graph.n_edges()); 
		INFO("GPB: K=%d, alpha=%.3lf, beta=%.3lf, %s\n", K, alpha, beta, (sample_deep ? "sample multinomial" : "skip multinomial sampling"));
		INFO("GPB: Fpshp=%.3lf, Fprte=%.3lf, Bpshp=%.3lf, Bprte=%.3lf\n", Fpshp, Fprte, Bpshp, Bprte);
		init();
	}

	void output_network_and_block();

    void gibbs(int burnin, int Ns);

    void save(const string& prefix, const mat& F, const mat& B) const;

	void load(const string& prefix);

	mat link_component();

	vector<set<string>> get_community(const mat& component, bool overlap=false);

	// return score, label pairs
    vector<pair<float,int>> link_prediction(const Graph::Heldout& test, bool save_in_file=true, const string& filename="heldout.txt") const;

	static int ECHO_PER_ITERS;

	static int SAVE_PER_ITERS;

	void set_dir(string dir) {this->save_dir = dir;}

private:
	Graph& graph;

	int N;

    int K;

    double Fpshp;

    double Fprte;

    double Bpshp;

    double Bprte;

	double beta;

	double alpha;

    mat Frte;

    mat Fshp;
    
    mat F;

    mat Bshp;

    mat Brte;

    mat B;

	mat phi;

	string save_dir;

	bool sample_deep;

    void init();

    double compute_elbo(const mat& F, const mat& B) const;

	//double heldout_likelihood(const Graph::Heldout& heldout) const;

    void sample_phi(int i, int j, bool sample_deep, bool accept_zero);

    void sample_F(const mat& Fshp, const mat& Frte);

    void sample_B(const mat& Bshp, const mat& Brte);

    default_random_engine generator;
};

#endif
