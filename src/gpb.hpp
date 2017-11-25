#ifndef GPB_H
#define GPB_H

#include <vector>
#include <armadillo>
#include <random>
#include "asa103.hpp"
#include "graph.hpp"
#include "utils.hpp"

using namespace std;
using namespace arma;

class GPB {

public:
    GPB(Graph& g, int k, double Fpshp=0.3, double Fprte=1.0, double Bpshp=0.3, double Bprte=1.0, const char* dir=".", bool sample_deep=false, int alpha=1):
    	graph(g), K(k), Fpshp(Fpshp), Fprte(Fprte), Bpshp(Bpshp), Bprte(Bprte), dir(dir), sample_deep(sample_deep) 
	{ 
		INFO("GPB settings: %s graph with %d nodes %d edges.\n", (graph.is_directed() ? "Directed" : "Undirected"), graph.n_nodes(), graph.n_edges()); 
		INFO("GPB settings: K(%d), alpha(%d), %s.\n", K, alpha, (sample_deep ? "sample deep" : "not sample deep"));
		INFO("GPB settings: Fpshp(%lf), Fprte(%lf), Bpshp(%lf), Bprte(%lf).\n", K, alpha, Fpshp, Fprte, Bpshp, Bprte);
		init(alpha);
	}

    void gibbs(int burnin, int Ns);

    // void vi(int max_iters);

    void save(const string& prefix, const mat& F, const mat& B) const;

	void load(const string& prefix);

	mat link_component();

	vector<set<string>> get_community(const mat& component, bool overlap=false);

	// if grid is true, do grid search
	// return best accuracy and set thresh to best thresh
	// if grid is false, return accuracy with thresh
    double link_prediction(const Graph::Heldout& test, bool grid, double& thresh) const;

	static int ECHO_PER_ITERS;

	static int SAVE_PER_ITERS;

	void set_dir(string dir) {this->dir = dir;}

private:
	const Graph& graph;

	int N;

    int K;

    double Fpshp;

    double Fprte;

    double Bpshp;

    double Bprte;

    mat Frte;

    mat Fshp;
    
    mat ElnF;

    mat F;

    mat ElnB;

    mat Bshp;

    mat Brte;

    mat B;

	string dir;

	bool sample_deep;

    void init(int alpha=1);

    // mat estimate_phi(int i, int j);

    double compute_elbo(const mat& F, const mat& B) const;


	//double heldout_likelihood(const Graph::Heldout& heldout) const;

    // void compute_exp_ln();

    mat sample_phi(int i, int j, bool sample_deep);

    void sample_F(const mat& Fshp, const mat& Frte);

    void sample_B(const mat& Bshp, const mat& Brte);

    default_random_engine generator;
};

#endif