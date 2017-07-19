#include "ppl.hpp"
#include "metrics.hpp"

void
PPL::init()
{
    if (this->initialized == true)
        return;

    generator.seed(rd());
    uniform_real_distribution<> dis(0.0, 1.0);

    // resize and resize init and init

    // random gamma
    G.resize(N);
    for (int i=0; i<N; ++i)
    {
        G[i].resize(topicnum);
        for (int k=0; k<topicnum; ++k)
            G[i][k] = dis(generator);
        normalize(G[i]);
    }
    // random W
    W.resize(topicnum);
    for (int i=0; i<topicnum; ++i)
    {
        W[i].resize(topicnum);
        for (int k=0; k<topicnum; ++k)
            W[i][k] = dis(generator);
        if (assortative)
            W[i][i] += 0.5;
    }
    normalize(W);
    // set A to outdegree
    A.resize(N);
    for (int i=0; i<N; ++i)
        A[i] = network[i].size();
    normalize(A);
    // set B to indegree
    B.resize(N, 0.0);
    for (int i=0; i<N; ++i)
        for (auto it=network[i].cbegin(); it!=network[i].cend(); ++it)
            B[*it]++;
    normalize(B);
    // resize Wt
    Wt.resize(topicnum);
    for (int i=0; i<topicnum; ++i)
        Wt[i].resize(topicnum);
    //resize Ni No
    Ni.resize(N);
    for (int i=0; i<N; ++i)
        Ni[i].resize(topicnum);
    No.resize(N);
    for (int i=0; i<N; ++i)
        No[i].resize(topicnum);

    // resize ETA, TAU
    ETA.resize(topicnum);
    TAU.resize(topicnum);
    update_eta_tau(); 

    converged.resize(N, -1);

    this->initialized = true;
}

void
PPL::run(int max_iters)
{
    init();

    vector<vector<double>> Q(topicnum, vector<double>(topicnum));

    int iter = 1;
    while (iter <= max_iters)
    {
        zero(Ni);
        zero(No);
        zero(Wt);

        for (int i=0; i<N; ++i)
            for (auto it=network[i].cbegin(); it!=network[i].cend(); ++it)
                E_step(Q, i, *it);

        M_step();
        update_eta_tau();
        cout << setfill('0') << setw(3) << iter << " " << likelihood() << flush << '\r';
        iter++;
    }
    cout << "final: " << likelihood() << endl;
}

void
PPL::update_eta_tau()
{
    zero(ETA);
    zero(TAU);
    for (int i=0; i<N; ++i)
        for (int k=0; k<topicnum; ++k)
        {
            ETA[k] += G[i][k]*A[i];
            TAU[k] += G[i][k]*B[i];
        }
}


void
PPL::E_step(vector<vector<double>>& Q, int i, int j)
{
    zero(Q);
    int p = converged[i];
    int q = converged[j];
    if (p != -1 && q != -1)
    {
        Q[p][q] = 1.0; // 都已经收敛了:)
    }
    else if (p != -1 && q == -1) // 收敛了一个:|
    {
        for (int h=0; h<topicnum; ++h)
            Q[p][h] = (G[i][p]*A[i]/ETA[p]) * (G[j][h]*B[j]/TAU[h]) * W[p][h];
        normalize(Q[p]);
    }
    else if (p == -1 && q != -1)
    {
        for (int g=0; g<topicnum; ++g)
            Q[g][q] = (G[i][g]*A[i]/ETA[g]) * (G[j][q]*B[j]/TAU[q]) * W[g][q];
        normalize(Q);
    }
    else if (q == -1 && p == -1) // 都没有收敛:(
    {
        for (int g=0; g<topicnum; ++g)
            for (int h=0; h<topicnum; ++h)
                Q[g][h] = (G[i][g]*A[i]/ETA[g]) * (G[j][h]*B[j]/TAU[h]) * W[g][h];
        normalize(Q);
    }
    else
    {
        cerr << "E_step exit abnormally!!";
        exit(1);
    }
    // update sufficient statistics
    for (int g=0; g<topicnum; ++g)
    {
        for (int h=0; h<topicnum; ++h)
        {
            No[i][g] += Q[g][h];
            Ni[j][h] += Q[g][h];
            Wt[g][h] += Q[g][h];
        }
    }
}

void
PPL::M_step()
{
    vector<double> ni(N, .0);
    vector<double> no(N, .0);
    vector<double> eta(topicnum, .0);
    vector<double> tau(topicnum, .0);
    for (int i=0; i<N; ++i)
        for (int j=0; j<topicnum; ++j)
        {
            ni[i] += Ni[i][j];
            no[i] += No[i][j];
        }
    for (int g=0; g<topicnum; ++g)
        for (int h=0; h<topicnum; ++h)
        {
            eta[g] += Wt[g][h];
            tau[h] += Wt[g][h];
        }
    for (int g=0; g<topicnum; ++g)
    {
        eta[g] /= ETA[g];
        tau[g] /= TAU[g];
    }

    // maximize G
    for (int i=0; i<N; ++i)
    {
        if (converged[i] != -1)
            continue;
        for (int k=0; k<topicnum; ++k)
            G[i][k] = (Ni[i][k] + No[i][k]) / (eta[k]*A[i]+tau[k]*B[i]);

        normalize(G[i]);

        for (int k=0; k<topicnum; ++k)
            if (G[i][k] >= 1-1e-4)
                converged[i] = k;
    }

    for (int i=0; i<N; ++i)
    {
        double sum = .0;
        for (int k=0; k<topicnum; ++k)
            sum += eta[k]*G[i][k];
        A[i] = no[i] / sum;
    }
    normalize(A);
    for (int i=0; i<N; ++i)
    {
        double sum = .0;
        for (int k=0; k<topicnum; ++k)
            sum += tau[k]*G[i][k];
        B[i] = ni[i] / sum;
    }
    normalize(B);

    normalize(Wt);
    swap(W, Wt);
}

double
PPL::likelihood() const
{
    double likelihood = .0;
    for (int i=0; i<N; ++i)
        for (auto it=network[i].cbegin(); it!=network[i].cend(); ++it)
        {
            double edge_likelihood = .0;
            for (int g=0; g<topicnum; ++g)
                for (int h=0; h<topicnum; ++h)
                    edge_likelihood += (G[i][g]*A[i]/ETA[g]) * (G[*it][h]*B[*it]/TAU[h]) * W[g][h];
            likelihood += log(edge_likelihood);
        }
    return likelihood;
}

void
PPL::save_model(string save_path) const
{
    string filename = save_path + "/w.txt";
    ofstream ofs(filename);
    for (int i=0; i<topicnum; ++i)
    {
        for (int j=0; j<topicnum; ++j)
            ofs << setprecision(4) << fixed << W[i][j] << " ";
        ofs << endl;
    }
    ofs.close();
    filename = save_path + "/g.txt";
    ofs.open(filename);
    for (int i=0; i<N; ++i)
    {
        ofs << id2str[i] << ": ";
        for (int j=0; j<topicnum; ++j)
            ofs << setprecision(4) << fixed << G[i][j] << " ";
        ofs << endl;
    }
    ofs.close();
}

vector<set<int>>
PPL::link_community(bool overlap)
{
    vector<set<int>> ss(topicnum);
    vector<vector<int>> cc(N, vector<int>(topicnum, 0));
    vector<vector<double>> q(topicnum, vector<double>(topicnum));
    double max; int maxi;
    for (int i=0; i<N; ++i)
        for (auto it=network[i].cbegin(); it!=network[i].cend(); ++it)
        {
            int j = *it;
            E_step(q, i, j);
            pair<int,int> pr;
            double max = .0;
            for (int k=0; k<topicnum; ++k)
            {
                auto largest = max_element(q[k].begin(), q[k].end());
                if (max < *largest)
                {
                    max = *largest;
                    pr = make_pair(k, largest-q[k].begin());
                }
            }
            cc[i][pr.first] += 1;
            cc[j][pr.second] += 1;
        }

    for (int i=0; i<N; ++i)
    {
        if (overlap)
        {
            for (int k=0; k<topicnum; ++k)
                if (cc[i][k] > 0)
                    ss[k].insert(i);
        }
        else
        {
            auto it = max_element(cc[i].begin(), cc[i].end());
            int k = it - cc[i].begin();
            ss[k].insert(i);
        }
    }
    return ss;
}

void usage()
{
    cout << "Usage: ppl filename #clusters #iters\n";
    cout << "       -g ground truth filename\n";
    cout << "       -u undirected graph\n";
    cout << "       -a assortative\n";

    exit(1);
}

void
PPL::normalize(vector<double>& v)
{
    double s = .0;
    int len = v.size();
    for (int i=0; i<len; ++i)
        s += v[i];
    for (int i=0; i<len; ++i)
        v[i] /= s;
}

void
PPL::normalize(vector<vector<double>>& v)
{
    double s = .0;
    int len = v.size();
    for (int i=0; i<len; ++i)
    {
        int len2 = v[i].size();
        for (int j=0; j<len2; ++j)
            s += v[i][j];
    }
    for (int i=0; i<len; ++i)
    {
        int len2 = v[i].size();
        for (int j=0; j<len2; ++j)
        {
            v[i][j] /= s;
        } 
    }
}

void
PPL::zero(vector<vector<double>>& v)
{
    int len = v.size();
    for (int i=0; i<len; ++i)
    {
        int len2 = v[i].size();
        for (int j=0; j<len2; ++j)
            v[i][j] = .0;
    }
}

void
PPL::zero(vector<double>& v)
{
    int len = v.size();
    for (int i=0; i<len; ++i)
        v[i] = .0;
}


int main(int argc, char* argv[])
{
    if (argc < 4)
        usage();

    bool directed = true;
    bool assort = false;
    char *gn = NULL;

    char *input = argv[1];

    int K = atoi(argv[2]);

    int max_iters = atoi(argv[3]);

    bool overlap = true;

    for (int i=4; i<argc; i++)
    {
        if (strcmp(argv[i], "-g")==0)
            gn = argv[++i];
        else if (strcmp(argv[i], "-u")==0)
            directed = false;
        else if (strcmp(argv[i], "-a")==0)
            assort = true;
        else if (strcmp(argv[i], "-disjoint")==0)
            overlap = false;
        else
            usage();
    }

    PPL *ppl = new PPL(K, assort);
    ppl->read_from_file(input, gn, directed);
    ppl->run(max_iters);

    auto g1 = ppl->ground_truth;
    auto g2 = ppl->link_community(overlap);
    ppl->save_community(g2);
    if (gn)
    {
        overlap ? Metrics<int>::ONMI(g1, g2) : Metrics<int>::NNMI(g1, g2);
        Metrics<int>::F1(g1,g2);
    }
    return 0;
}
