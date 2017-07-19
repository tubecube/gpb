#ifndef CODA_H
#define CODA_H

#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <random>
#include <armadillo>
#include <iomanip>
#include <iostream>
#include "graph.hpp"

using namespace std;
using namespace arma;

class CODA : public Graph {

typedef double (CODA::*likelihoodFun) (const rowvec&, int) const;

typedef rowvec (CODA::*gradientFun) (int) const;

public:

    CODA(int k):K(k) {}

    void run(int max_iters);

    vector<set<int>> link_community(bool overlap, int type);

    vector<set<int>> node_community(bool overlap, int type);

private:

    const int K;

    mat F; // mat = Mat<double>

    mat H;

    rowvec Hsum; // rowvec = Row<double>

    rowvec Fsum;

    rowvec Fgradient(int) const;

    rowvec Hgradient(int) const;

    double Flikelihood(const rowvec&, int) const;

    double Hlikelihood(const rowvec&, int) const;

    void BLS(mat&, rowvec&, int, gradientFun, likelihoodFun);

    void init();

    double likelihood() const;
};

#endif
