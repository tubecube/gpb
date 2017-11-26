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
		// pair unhashable
		typedef set<pair<int,int>> pair_set;
		pair_set pairs[2];
	};

	vector<Heldout> heldouts;

	const Heldout* push_heldout(int N0, int N1);

	const Heldout* push_heldout_with_percentage(double pctg0, double pctg1)
	{
		if (pctg0 < 1 && pctg0 > 0 && pctg1 < 1 && pctg1 > 0)
			return push_heldout(int(Nzeros*pctg0), int(Nones*pctg1));
		else
		{
			ERROR("Heldout: not a valid percentage!\n");
			return NULL;
		}
	}

	int pop_heldout();

	bool check_in_heldouts(int source, int dest, bool link) const;

	void check_edge(int& source, int& dest) const
	{
		if (!directed && source > dest)
			std::swap(source, dest);
	}

	int N;
	unsigned long Nones;
	unsigned long Nzeros;

private:
	double ratio0;
	double ratio1;
	bool directed;
	vector<set<int>> edges;
	unordered_map<string, int> str2id;
	vector<string> id2str;
    int check_id(const string&);
};

#endif
