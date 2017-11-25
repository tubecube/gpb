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
	cout << "-grid: use grid search for link prediction\n";
	cout << "-deep: sample deep\n";
	cout << "-dir ?: directory to save model(default: model_k)\n";
	cout << "-alpha ?: factor for block matirx's diagnal shape priors(defalut: 1)\n";
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
	bool sample_deep = false;
	bool grid = false;
	string ground_truth="";
	Metrics<string>::file_type type = Metrics<string>::type1;
	int K = 10;
    int Ns = 50;
    int burnin = 50;
	int Alpha = 1;
    double Fpshp = 0.3;
    double Fprte = 1.0;
    double Bpshp = 0.3;
    double Bprte = 1.0;
    string save_dir="";
	Graph::Heldout* hp = NULL;

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

		else if (strcmp(argv[i], "-grid")==0)
			grid = true;

		else if (strcmp(argv[i], "-deep")==0)
			sample_deep = true;

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
			Alpha = atoi(argv[++i]);

		else {
			usage();
			exit(1);
		}
    }

	if (save_dir.size() == 0)
		save_dir = "model_" + to_string(K);
	makedir(save_dir.c_str());

	Graph graph = Graph(input, directed);

	if (heldout)
	{
		int n0[1] = {100};
		int n1[1] = {100};
		int *sizes[2] = {n0, n1};
		hp = graph.create_heldouts(sizes, 1);
	}

    GPB gpb(graph, K, Fpshp, Fprte, Bpshp, Bprte, save_dir.c_str(), sample_deep, Alpha);

    if (Gibbs)	gpb.gibbs(burnin, Ns);
    // else	gpb->vi(vimax);

	if (heldout)
	{
		double thresh = hp[0].ratio1;
		gpb.link_prediction(hp[0], grid, thresh);
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
