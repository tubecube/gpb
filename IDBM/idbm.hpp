#include <utility>
#include <vector>
#include <random>
#include <iomanip>
#include <cstdlib>
#include "Graph.h"

using namespace std;

class IDBM : public Graph {

public:
    IDBM(int k, double a, double b):topicnum(k), alpha(a), beta(b), initialized(false) {}

    void gibbs(int burnin); 

    vector<set<string>> link_community() const;
    vector<set<string>> community() const;
    void save_PSI(string file);

private:

    int topicnum;
    double alpha;
    double beta;
    double baselh;
    bool initialized;
    int iter;

    vector<vector<int>> N;
    vector<vector<int>> P;
    vector<vector<int>> R;
    vector<int> Nt;
    vector<int> Pt;
    vector<pair<int,int>> Z;
    // 随机数发生器
    default_random_engine generator;
    // 产生随机数发生器的种子
    random_device rd;

    void init();

    pair<int,int> sampling(pair<int,int> nodepr, pair<int,int> oldpr);
    pair<int,int> sampling2(pair<int,int> nodepr, pair<int,int> oldpr);

    double likelihood() const;

};
