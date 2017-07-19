#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <set>
#include <map>
#include <ctime>

using namespace std;

class Graph {

public:
    Graph():N(0), M(0){}

    Graph(const string &filename, const string &ground, bool directed) {
        read_from_file(filename, ground, directed);
    }

    void read_from_file(const string &filename, const string &ground, bool directed=true, bool heldout=false);

    int N;//nodenum

    int M;//edgenum
    
    vector<string> id2str;

    unordered_map<string, int> str2id;

    bool directed;

    vector<set<int>> network;

    vector<set<int>> network2;

    map<pair<int,int>, bool> heldout;

    vector<set<int>> ground_truth;

    set<int> ground_set;

    void add_edge(string start, string end);
    
    int get_id(const string& str);

    map<pair<int,int>, bool>& create_heldout(size_t size1, size_t size2);

    virtual void run(int max_iters) = 0;

    void cut(vector<set<int>>&);

    void save_community(vector<set<int>>&);

    void create_ground_truth(const string&);
};

#endif
