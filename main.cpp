#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include "metrics.hpp"
#include "gpb.hpp"

using namespace std;

void usage()
{
    cout << "gpb    \"input_edge_file\" \"#clusters\" \"ground_truth_file(optional)\"\n";
    cout << "other optional parameters:\n";
    cout << "       -undirected ------------ undirected graph\n";
    cout << "       -heldout -------------- use heldout set \n";
    cout << "       -disjoint ------------- create disjoint communities\n";
    cout << "       -thresh --------------- thresh for detecting overlapping comms\n";
    cout << "       -Fpshp,Fprte,Bpshp,Bprte ------- hyper-parameters\n";
    cout << "       -vi ------------------- variational inference\n";
    cout << "       -max N -------- variational inference max iterations\n";
    cout << "       -burnin N ----------- mcmc burnin\n";
    cout << "       -sample N ----------- mcmc samples\n";
    exit(1);
}

int main(int argc, char* argv[]) {
    if (argc < 3)
        usage();

    string input = argv[1];
    int K = atoi(argv[2]);

    string ground = "";

    if (argc > 3) {
        string argv3 = argv[3];
        if (argv3[0] != '-')
            ground = argv3;
    }

    bool directed = true;
    bool heldout = false;
    int max_iters = -1;
    double thresh = 0;
    bool overlap = true;
    double Fpshp = 0.3;
    double Fprte = 1.0;
    double Bpshp = 0.3;
    double Bprte = 1.0;
    string method = "mcmc";
    string save_dir;
    int sample = 50;
    int burnin = 100;

    for (int i=3; i<argc; ++i)
    {
        if (strcmp(argv[i], "-max")==0)
            max_iters = atoi(argv[++i]);
        else if (strcmp(argv[i], "-thresh")==0)
            thresh = atof(argv[++i]);
        else if (strcmp(argv[i], "-disjoint")==0)
            overlap = false;
        else if (strcmp(argv[i], "-undirected")==0)
            directed = false;
        else if (strcmp(argv[i], "-heldout")==0)
            heldout = true;
        else if (strcmp(argv[i], "-Fpshp")==0)
            Fpshp = atof(argv[++i]);
        else if (strcmp(argv[i], "-Fprte")==0)
            Fprte = atof(argv[++i]);
        else if (strcmp(argv[i], "-Bpshp")==0)
            Bpshp = atof(argv[++i]);
        else if (strcmp(argv[i], "-Bprte")==0)
            Bprte = atof(argv[++i]);
        else if (strcmp(argv[i], "-vi")==0)
            method = "vi";
        else if (strcmp(argv[i], "-burnin")==0)
            burnin = atoi(argv[++i]);
        else if (strcmp(argv[i], "-sample")==0)
            sample = atoi(argv[++i]);
        else if (strcmp(argv[i], "-save")==0)
            save_dir = argv[++i];
    }

    GPB *model = new GPB(K, Fpshp, Fprte, Bpshp, Bprte);
    model->read_from_file(input, ground, directed, heldout);
    if (method == "mcmc")
        model->mcmc(burnin, sample);
    else if (method == "vi")
        model->run(max_iters);
    // bool overlap = Metrics<int>::is_overlap(model->ground_truth);
    if (save_dir.size())
        model->save_model(save_dir);
    auto g1 = model->ground_truth;
    auto g2 = model->node_community(overlap, thresh);
    model->save_community(g2);
    /*
    if (g1.size() > 0) {
        model->cut(g2);
        overlap ? Metrics<int>::ONMI(g1,g2) : Metrics<int>::NNMI(g1, g2);
        Metrics<int>::F1(g1, g2);
    }
    */
    return 0;
}
