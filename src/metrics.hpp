#ifndef MEstringRICS_H
#define MEstringRICS_H

#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <algorithm>
#include "utils.hpp"

using namespace std;

class Metrics {
public:
	/*
	type1: each line indicates a community, with each node seperated by tab.
	type2: each line indicates a node community pair, seperated by tab.
	*/
	typedef enum {type1, type2} file_type;
	/* from community file to community set */
	static vector<set<string>> file_to_set(const string& filename, file_type type=type1)
	{
		if (type == type1)
			return file_to_set1(filename);
		else
			return file_to_set2(filename);
	}
	/* from community set to community file */
	static void set_to_file(const string& filename, const vector<set<string>>& community, file_type type=type1)
	{
		if (type == type1)
			set_to_file1(filename, community);
		else
			set_to_file2(filename, community);
	}
	/* NMI score of disjoint communities */
	static double NNMI(const vector< set<string> >& g1, const vector< set<string> >& g2);
	/* NMI score of overlapping communities
	   MacDaid(2011). Normalized mutual information to evaluate overlapping community finding algorithms.*/	
	static vector<double> ONMI(const vector< set<string> >& g1, const vector< set<string> >& g2);
	/* F1 score */
	static double F1(const vector< set<string> >& g1, const vector< set<string> >& g2);
	/* check if communities overlap */
	static bool is_overlap(const vector< set<string> >& g);
	/* keep only overlapping nodes in two communitiy sets */
	static void keepOnlyOverlap(vector< set<string> >& g1, vector< set<string> >& g2);
	static double F12(const vector<set<string>>& g1, const vector<set<string>>& g2);
private:
	static vector<set<string>> file_to_set1(const string& filename);
	static vector<set<string>> file_to_set2(const string& filename);
	static void set_to_file1(const string& filename, const vector<set<string>>& community);
	static void set_to_file2(const string& filename, const vector<set<string>>& community);
	static size_t overlap_size(const set<string>&, const set<string>&);
	static void onmi_one_turn(const vector< set<string> >&, const vector< set<string> >&, size_t, double&, double&, double&);
	static double F1_one_turn(const vector< set<string> >&, const vector< set<string> >&);
	static double H(double p);
};

#endif
