#ifndef HPF_H
#define HPF_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <list>
#include <set>
#include "asa103.hpp"
#include "graph.hpp"

using namespace arma;
using namespace std;

class HPF : public Graph {

public:
    HPF(int k, double __a=0.3, double __b=0.3, double __c=0.3, double __d=0.3, double _a=0.3, double _c=0.3):K(k),__a(__a),__b(__b),__c(__c),__d(__d),_a(_a),_c(_c){}
    void run(int max_iters);
    vector<set<int>> link_community(bool overlap, double thresh);
    vector<set<int>> node_community(bool overlap, int type);
    void save_model(string dirname) const;
    void load_model(string dirname);

private:
    int K;
    // hyper
    double __a; //0.3
    double __b; // def:0.3
    double __c; //0.3
    double __d; // def:0.3
    double _a; //0.3
    double _c; //0.3
    
    mat lrte; // lambda rate (K*N)
    mat lshp; // lambda shape (K*N)
    mat grte;
    mat gshp;

    rowvec krte; // kappa rte
    rowvec kshp; // kappa shape
    rowvec trte; // tau rte
    rowvec tshp; // tau shp

    mat lnrte;
    mat lnshp;
    mat gnrte;
    mat gnshp;

    mat Elngamma;
    mat Elnlambda;

    vector<list<int>> gactive;
    vector<list<int>> lactive;

    uvec gconv;
    uvec lconv;

    void init();

    vec compute_phi(int u, int i);

    double likelihood() const;

    double validation_likelihood() const;

    void compute_exp_ln();

    void check_converge();
};

#endif
