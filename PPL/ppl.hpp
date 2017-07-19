#ifndef PPL_H
#define PPL_H

#include <vector>
#include <iomanip>
#include <random>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "graph.hpp"

using namespace std;

class PPL : public Graph {

public:
    PPL(int k, bool a):topicnum(k),assortative(a),initialized(false){}

    void run(int max_iters);

    void save_model(string path) const;

    vector<set<int>> link_community(bool overlap);

private:

    int topicnum;
    int max_iters;
    bool initialized;
    bool assortative;

    vector<vector<double>> G;
    vector<vector<double>> W;
    vector<vector<double>> Ni;
    vector<vector<double>> No;
    vector<vector<double>> Wt;
    /*
       A is productivity;
       B is pupularity;
    */
    vector<double> A;
    vector<double> B;
    vector<double> ETA;
    vector<double> TAU;
    vector<int> converged;

    void init();

    void update_eta_tau();

    void E_step(vector<vector<double>>& q, int i, int j);

    void M_step();

    double likelihood() const;

    void normalize(vector<double>& v);

    void normalize(vector<vector<double>>& v);

    void zero(vector<vector<double>>& v);

    void zero(vector<double>& v);

    default_random_engine generator;

    random_device rd;
};

#endif
