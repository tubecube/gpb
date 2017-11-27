#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include "roc.hpp"
#include "utils.hpp"
#include "graph.hpp"
#include "metrics.hpp"
#include "gpb.hpp"

using namespace std;

void usage();

void compare_with_ground_truth(vector<set<string>>& community, vector<set<string>>& ground_truth);

int main(int argc, char* argv[]) {
	Logger::_LEVEL = Logger::INFO;
    if (argc <= 1) {
        usage();
		exit(1);
	}
	// 输入文件
    string input = argv[1];
	// default settings
    bool directed = true;
	bool sample_deep = true;
	bool lp = false;
	bool comp_grnd = false;
	string grnd_file = "";
	Metrics::file_type type = Metrics::type1;
	int K = 10;
    int Ns = 50;
    int burnin = 50;
	double alpha = 1.0;
	double beta = 0.0;
    double Fpshp = 0.3;
    double Fprte = 1.0;
    double Bpshp = 0.3;
    double Bprte = 1.0;
    string save_dir="";

    for (int i=2; i<argc; ++i)
    {
		// number of dimensions
		if (strcmp(argv[i], "-k")==0)
			K = atoi(argv[++i]);
		// burnin in Gibbs sampling
        else if (strcmp(argv[i], "-burn")==0)
            burnin = atoi(argv[++i]);
		// #samples in Gibbs sampling
        else if (strcmp(argv[i], "-sample")==0)
            Ns = atoi(argv[++i]);
		// indicate an undirected graph
        else if (strcmp(argv[i], "-u")==0)
            directed = false;
		// do link prediction
        else if (strcmp(argv[i], "-lp")==0)
            lp = true;
		// not sample_deep
		else if (strcmp(argv[i], "-nsd")==0)
			sample_deep = false;
		// hyperparameters
        else if (strcmp(argv[i], "-Fpshp")==0)
            Fpshp = atof(argv[++i]);
        else if (strcmp(argv[i], "-Fprte")==0)
            Fprte = atof(argv[++i]);
        else if (strcmp(argv[i], "-Bpshp")==0)
            Bpshp = atof(argv[++i]);
        else if (strcmp(argv[i], "-Bprte")==0)
            Bprte = atof(argv[++i]);
		// saved directory
        else if (strcmp(argv[i], "-dir")==0)
            save_dir = argv[++i];
		// compare with ground truth
		else if (strcmp(argv[i], "-grnd")==0)
		{
			grnd_file = argv[++i];
			comp_grnd = true;
		}
		// ground truth file type
		else if (strcmp(argv[i], "-type2")==0)
			type = Metrics::type2;
		//
		else if (strcmp(argv[i], "-alpha")==0)
			alpha = atof(argv[++i]);
		//
		else if (strcmp(argv[i], "-beta")==0)
			beta = atof(argv[++i]);
		else {
			usage();
			exit(1);
		}
    }

	if (save_dir.size() == 0)
		save_dir = "gpb_K" + to_string(K);

	Graph graph = Graph(input, directed);

	const Graph::Heldout* test;

	if (lp)
	{
		double ratio = 0.1;
		test = graph.push_heldout_with_percentage(ratio, ratio);
	}
    GPB gpb(graph,K,Fpshp,Fprte,Bpshp,Bprte,save_dir.c_str(),sample_deep,alpha,beta);
    gpb.gibbs(burnin, Ns);
	if (lp)
	{
		vector<pair<float,int>> data = gpb.link_prediction(*test);
		ROC roc(data);
		INFO("Link prediction AUC: %f\n", roc.getAreaUnderCurve());
		roc.writeToFile(save_dir+"/roc-"+input);
	}
	vector<set<string>> community = gpb.get_community(gpb.link_component());
	Metrics::set_to_file(save_dir+"/community.dat", community, type);
	if (comp_grnd)
	{
		vector<set<string>> ground_truth = Metrics::file_to_set(grnd_file, type);
		compare_with_ground_truth(community, ground_truth);
		Metrics::set_to_file(save_dir+"/ground_truth.dat", ground_truth, type);
	}
    return 0;
}

void usage()
{
    cout << "USAGE: gpb INPUT_FILE\n";
	cout << "\n";
    cout << "OPTIONAL parameters:\n";
    cout << "-u: input graph is undirected\n";
	cout << "-k: set feature vector dimension(default: 10)\n";
	cout << "-grnd: set comm filename\n";
	cout << "-type2: each line is a node-comm pair(default: each line is a comm)\n";
    cout << "-burn: set gibbs sampling burnin(default: 50)\n";
    cout << "-sample: set num of samples(default: 50)\n";
    cout << "-lp: do link prediction\n";
	cout << "-dir: set directory to save\n";
	cout << "-alpha: set factor for homophily defalut: 1.0\n";
	cout << "-beta: set additional samples per iter default: 0.0\n";
	cout << "-nsd: not sample deep\n";
    cout << "-Fpshp: set hyper parameters default: 0.3\n";
    cout << "-Fprte: set hyper parameters default: 1.0\n";
    cout << "-Bpshp: set hyper parameters default: 0.3\n";
    cout << "-Bprte: set hyper parameters default: 1.0\n";
}

void compare_with_ground_truth(vector<set<string>>& community, vector<set<string>>& ground_truth)
{
	Metrics::keepOnlyOverlap(community, ground_truth);
	INFO("Community detection metrics:\n");
	Metrics::F1(community, ground_truth); 
	Metrics::NNMI(community, ground_truth);
	Metrics::ONMI(community, ground_truth);
}
