#ifndef METRICS_H
#define METRICS_H

#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>
#include "utils.hpp"

using namespace std;

template <class T>
class Metrics {
public:
	typedef enum {type1, type2} file_type;
	static vector<set<T>> file_to_set(const string& filename, file_type type=type1)
	{
		if (type == type1)
			return file_to_set1(filename);
		else
			return file_to_set2(filename);
	}
	static void set_to_file(const string& filename, const vector<set<T>>& community, file_type type=type1)
	{
		if (type == type1)
			set_to_file1(filename, community);
		else
			set_to_file2(filename, community);
	}
    static double NNMI(const vector< set<T> >& g1, const vector< set<T> >& g2);
    static vector<double> ONMI(const vector< set<T> >& g1, const vector< set<T> >& g2);
    static double F1(const vector< set<T> >& g1, const vector< set<T> >& g2);
    static bool is_overlap(const vector< set<T> >& g);
    static void keepOnlyOverlap(vector< set<T> >& g1, vector< set<T> >& g2);
private:
	// each row is a community
	static vector<set<T>> file_to_set1(const string& filename);
	// each row is a node-community pair
	static vector<set<T>> file_to_set2(const string& filename);
	static void set_to_file1(const string& filename, const vector<set<T>>& community);
	static void set_to_file2(const string& filename, const vector<set<T>>& community);
    static size_t overlap_size(const set<T>&, const set<T>&);
    static void onmi_one_turn(const vector< set<T> >&, const vector< set<T> >&, size_t, double&, double&, double&);
    static double F1_one_turn(const vector< set<T> >&, const vector< set<T> >&);
    static double H(double p);
};

#include "metrics.cpp"

#endif
/*
void usage()
{
    cout << " metric file1 file2 [1(for F1)] [2(for ONMI)] [3(for NNMI)] [-c(keep overlap part)]" << endl;
    exit(1);
}

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        usage();
    }

    if (strcmp(argv[1], "check") == 0)
    {
        char *file = argv[2];
        cout << Metrics<string>::is_overlap(fileToSet(file)) << endl;
        return 0;
    }

    char *file1 = argv[1];
    char *file2 = argv[2];

    auto g1 = fileToSet(file1);
    auto g2 = fileToSet(file2);

    int method;

    if (argc > 3)
        method = atoi(argv[3]);
    else
        usage();

    if (argc > 4)
        if (strcmp(argv[4], "-c") == 0)
            Metrics<string>::keepOnlyOverlap(g1, g2);

    if (method == 1)
        Metrics<string>::F1(g1, g2);
    else if (method == 2)
        Metrics<string>::ONMI(g1, g2);
    else if (method == 3)
        Metrics<string>::NNMI(g1, g2);
    else
        usage();

    return 0;
}    
*/
