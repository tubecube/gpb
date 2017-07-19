#ifndef GPB_H
#define GPB_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <set>
#include <list>
#include <random>
#include "asa103.hpp"
#include "graph.hpp"

using namespace std;
using namespace arma;

class GPB : public Graph {

public:
    GPB(int k, double Fpshp=0.3, double Fprte=1.0, double Bpshp=0.3, double Bprte=1.0)
    :K(k), Fpshp(Fpshp), Fprte(Fprte), Bpshp(Bpshp), Bprte(Bprte){}

    void run(int max_iters);

    vector<set<int>> link_community(bool, double);

    vector<vector<int>> node_community(bool, double);

    void save_model(string dirname) const;

    void mcmc(int burnin, int sample);

    void save_community(vector<vector<int>>& comm);

private:
    int K;

    double Fpshp;

    double Fprte;

    double Bpshp;

    double Bprte;

    mat Frte; // 

    mat Fshp;
    
    mat ElnF;
    
    mat EF;

    mat F;

    mat EB; //

    mat ElnB;

    mat Bshp;

    mat Brte;

    mat B;

    void init();

    // sp_mat estimate_phi(int i, int j);

    mat estimate_phi(int i, int j);

    double likelihood() const;

    double validation_likelihood() const;

    void compute_exp_ln();

    sp_imat Network;

    mat sample_phi(int i, int j);

    mat sample_F(mat Fshp, mat Frte);

    mat sample_B(mat Bshp, mat Brte);

    int iter;

    default_random_engine generator;

    random_device rd;    
};

#endif
