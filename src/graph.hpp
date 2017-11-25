#ifndef GRAPH_H
#define GRAPH_H

#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <ctime>
#include <armadillo>
#include "utils.hpp"

using namespace std;
using namespace arma;

class Graph
{
public:
	SpMat<unsigned> network;

    Graph():N(0),Nones(0),Nzeros(0){}

    Graph(const string &filename, bool directed):N(0),Nones(0),Nzeros(0)
	{
        read_from_file(filename, directed);
    }

    int read_from_file(const string &filename, bool directed);

	bool is_directed() const {return directed;}

	int n_nodes() const {return N;}

	unsigned n_edges() const {return Nones;}

	string get_str(int id) const { return id2str[id]; }

	int get_id(const string& str) const
	{
		auto iter = str2id.find(str);
		if (iter == str2id.end())
			return -1;
		return iter->second;
	}

	struct Heldout
	{
		typedef set<pair<int,int>> pair_set;
		Heldout():ones(pairs[1]),zeros(pairs[0]) {}
		pair_set pairs[2];
		pair_set &ones;
		pair_set &zeros;
		double ratio0;
		double ratio1;
	};

	/* create N heldout sets: */
	/* sizes[0] is # of nonlinks in each heldout set. */
	/* sizes[1] is # of links in each heldout set. */
	Heldout* create_heldouts(int *sizes[2], int N);

private:
	int N;
	unsigned long Nones;
	unsigned long Nzeros;
	double ratio0;
	double ratio1;
	bool directed;
	vector<set<int>> edges;
	unordered_map<string, int> str2id;
	vector<string> id2str;
    int check_id(const string&);
};

#endif
