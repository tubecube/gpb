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

void usage()
{
    cout << "USAGE: gpb INPUT_FILE\n";
	cout << "\n";
    cout << "OPTIONAL parameters:\n";
	cout << "-k ?: feature vector dimensions(default: 10)\n";
	cout << "-ground ?: ground truth filename(default: None)\n";
	cout << "-type2: ground truth file with type2\n";
    cout << "-burn ?: Gibbs sampling burnin(default: 50)\n";
    cout << "-ns ?: samples in Gibbs sampling(default: 50)\n";
    cout << "-u: indicate an undirected graph\n";
    cout << "-test: use a heldout set for link prediction\n";
	cout << "-dir ?: directory to save model(default: model_k)\n";
	cout << "-alpha ?: factor for block matirx's diagnal shape priors(defalut: 1)\n";
	cout << "-beta ?: factor for additional sample zeros(default: 0.2)\n";
	cout << "-ndeep: not sample deep\n";
    cout << "-Fpshp ?: hyper parameters(default: 0.3)\n";
    cout << "-Fprte ?: hyper parameters(default: 1.0)\n";
    cout << "-Bpshp ?: hyper parameters(default: 0.3)\n";
    cout << "-Bprte ?: hyper parameters(default: 1.0)\n";
}


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
    bool heldout = false;
	bool sample_deep = true;
	string ground_truth="";
	Metrics<string>::file_type type = Metrics<string>::type1;
	int K = 10;
    int Ns = 50;
    int burnin = 50;
	double alpha = 1.0;
	double beta = 0.2;
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
        else if (strcmp(argv[i], "-ns")==0)
            Ns = atoi(argv[++i]);

		// indicate an undirected graph
        else if (strcmp(argv[i], "-u")==0)
            directed = false;

		// use a heldout set
        else if (strcmp(argv[i], "-test")==0)
            heldout = true;

		else if (strcmp(argv[i], "-ndeep")==0)
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

		else if (strcmp(argv[i], "-ground")==0)
			ground_truth = argv[++i];

		else if (strcmp(argv[i], "-type2")==0)
			type = Metrics<string>::type2;

		else if (strcmp(argv[i], "-alpha")==0)
			alpha = atof(argv[++i]);

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

	if (heldout)
	{
		int sign = graph.push_heldout(100, 1000);
		if (sign != 0)
			exit(sign);
	}

    GPB gpb(graph, K, Fpshp, Fprte, Bpshp, Bprte, save_dir.c_str(), sample_deep, alpha, beta);

    gpb.gibbs(burnin, Ns);

	if (heldout)
	{
		vector<pair<float,int>> data = gpb.link_prediction(graph.heldouts.back());
		ROC roc(data);
		INFO("AUC: %f\n", roc.getAreaUnderCurve());
	}

	vector<set<string>> community = gpb.get_community(gpb.link_component());
	Metrics<string>::set_to_file(save_dir+"/community.dat", community);
	if (ground_truth.size() > 0)
	{
		vector<set<string>> base = Metrics<string>::file_to_set(ground_truth, type);
    	Metrics<string>::F1(community, base);
		Metrics<string>::NNMI(community, base);
		Metrics<string>::ONMI(community, base);
	}
    return 0;
}
