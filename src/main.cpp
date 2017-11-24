#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include "utils.hpp"
#include "graph.hpp"
#include "metrics.hpp"
#include "gpb.hpp"

using namespace std;

void usage()
{
    cout << "USAGE: gpb INPUT_FILE [-u][-h][-v][-k ?][-m ?][-b ?][-s ?][-Fpshp ?][-Fprte ?][-Bpshp ?][-Bprte ?]\n";
	cout << "\n";
    cout << "OPTIONAL parameters:\n";
    cout << "-u: indicate an undirected graph\n";
    cout << "-h: use a heldout set\n";
    cout << "-v: use variational inference\n";
	cout << "-k ?: input dimensions(default: 10)\n";
    cout << "-m ?: max iterations in vi(default: 100)\n";
    cout << "-b ?: burnin in Gibbs sampling(default: 50)\n";
    cout << "-n ?: samples in Gibbs sampling(default: 50)\n";
	cout << "-s ?: directory to save model(default: model_k)\n";
	cout << "-g ?: ground truth filename(default: None)\n";
    cout << "-Fpshp ?: hyper parameters(default: 0.3)\n";
    cout << "-Fprte ?: hyper parameters(default: 1.0)\n";
    cout << "-Bpshp ?: hyper parameters(default: 0.3)\n";
    cout << "-Bprte ?: hyper parameters(default: 1.0)\n";
}

int main(int argc, char* argv[]) {
    if (argc <= 1) {
        usage();
		exit(1);
	}

	// 输入文件
    string input = argv[1];

	// default settings
    bool directed = true;
    bool heldout = false;
	bool Gibbs = true;
	int K = 10;
    int vimax = 100;
    int Ns = 50;
    int burnin = 50;
    double Fpshp = 0.3;
    double Fprte = 1.0;
    double Bpshp = 0.3;
    double Bprte = 1.0;
    string save_dir="";
	string ground_truth="";

    for (int i=2; i<argc; ++i)
    {
		// number of dimensions
		if (strcmp(argv[i], "-k")==0)
			K = atoi(argv[++i]);

		// use variational inference instead of Gibbs sampling
		else if (strcmp(argv[i], "-v")==0)
			Gibbs = false;

		// variational inference max iterations
		else if (strcmp(argv[i], "-m")==0)
            vimax = atoi(argv[++i]);

		// burnin in Gibbs sampling
        else if (strcmp(argv[i], "-b")==0)
            burnin = atoi(argv[++i]);

		// #samples in Gibbs sampling
        else if (strcmp(argv[i], "-n")==0)
            Ns = atoi(argv[++i]);

		// indicate an undirected graph
        else if (strcmp(argv[i], "-u")==0)
            directed = false;

		// use a heldout set
        else if (strcmp(argv[i], "-h")==0)
            heldout = true;

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
        else if (strcmp(argv[i], "-s")==0)
            save_dir = argv[++i];

		else if (strcmp(argv[i], "-g")==0)
			ground_truth = argv[++i];
		else {
			usage();
			exit(1);
		}
    }

	if (save_dir.size() == 0)
		save_dir = "model_" + to_string(K);
	makedir(save_dir.c_str());

	Graph graph = Graph(input, directed);

    GPB *gpb = new GPB(graph, K, Fpshp, Fprte, Bpshp, Bprte, save_dir.c_str());

    if (Gibbs)	gpb->gibbs(burnin, Ns);
    // else	gpb->vi(vimax);

	vector<set<string>> community = gpb->get_community(gpb->link_component());
	Metrics<string>::set_to_file(save_dir+"/community.dat", community);
	if (ground_truth.size() > 0)
	{
		vector<set<string>> base = Metrics<string>::file_to_set2(ground_truth);
    	Metrics<string>::F1(community, base);
	}
    return 0;
}
