#include "coda.h"

rowvec CODA::Fgradient(int u) const
{
    rowvec gradient = zeros<rowvec>(K);
    rowvec Hs(Hsum - H.row(u));
    
    for (auto itr = network[u].begin(); itr != network[u].end(); ++itr)
    {
        int v = *itr;
        double factor = dot(F.row(u),H.row(v));
        factor = exp(-factor) / (1.0 - exp(-factor));
        gradient += factor * H.row(v);
        Hs -= H.row(v);
    }
    gradient -= Hs;
    return gradient;
}

rowvec CODA::Hgradient(int v) const
{
    rowvec gradient = zeros<rowvec>(K);
    rowvec Fs(Fsum - F.row(v));
    
    for (auto itr = network2[v].begin(); itr != network2[v].end(); ++itr)
    {
        int u = *itr;
        double factor = dot(F.row(u),H.row(v));
        factor = exp(-factor) / (1.0 - exp(-factor));
        gradient += factor * F.row(u);
        Fs -= F.row(u);
    }
    gradient -= Fs;
    return gradient;
}

void CODA::run(int max_iters)
{
    init();
    double lh = .0;
    int iter = 1;
    for (; max_iters==-1 || iter<=max_iters; ++iter)
    {
        // update F
        for (int u=0; u<N; ++u)
            BLS(F, Fsum, u, &CODA::Fgradient, &CODA::Flikelihood);
        // update H
        for (int v=0; v<N; ++v)
            BLS(H, Hsum, v, &CODA::Hgradient, &CODA::Hlikelihood);
        double last = lh;
        lh = likelihood();
        cout << setfill('0') << setw(3) << iter << " " << lh << flush << '\r';
        if (abs( (lh-last)/lh ) < 1e-5)
            break;
    }
    cout << setfill('0') << setw(3) << iter << " " << lh << endl;
}

void CODA::BLS(mat& m, rowvec& sum, int u, gradientFun gfun, likelihoodFun lfun)
{
    rowvec grad = (this->*gfun)(u);
    grad = normalise(grad);
    double lh = (this->*lfun)(m.row(u), u);
    // double sqr = dot(grad,grad);
    double t = 1.0;
    double a = 0.5;
    int n = 3;
    rowvec newvec;
    do {
        newvec = m.row(u) + t*grad;
        newvec.elem( find(newvec < 0) ).zeros();
    } while (((this->*lfun)(newvec, u) - lh < a*t) && n-- && (t*=0.8));
    if (n == -1) return;
    sum -= m.row(u);
    m.row(u) = newvec;
    sum += newvec;
}

double CODA::Flikelihood(const rowvec& Fu, int u) const
{
    double likelihood = .0;

    rowvec Hs(Hsum - H.row(u));

    for (auto itr = network[u].begin(); itr != network[u].end(); ++itr)
    {
        int v = *itr;
        likelihood += log(1 - exp(-dot(Fu,H.row(v))));
        Hs -= H.row(v);
    }

    likelihood -= dot(Fu,Hs);
    return likelihood;
}

double CODA::Hlikelihood(const rowvec& Hv, int v) const
{
    double likelihood = .0;
    rowvec Fs(Fsum - F.row(v));

    for (auto itr = network2[v].begin(); itr != network2[v].end(); ++itr)
    {
        int u = *itr;
        likelihood += log(1 - exp(-dot(F.row(u), Hv)));
        Fs -= F.row(u);
    }

    likelihood -= dot(Hv,Fs);
    return likelihood;
}

void CODA::init()
{
    arma_rng::set_seed_random();
    F.randu(N,K);
    Fsum = sum(F,0);
    H.randu(N,K);
    Hsum = sum(H,0);
}

vector<set<int>>
CODA::link_community(bool overlap, int type)
{
    vector<set<int>> ss;
    return ss;
}

vector<set<int>>
CODA::node_community(bool overlap, int type)
{
    vector<set<int>> ss(K);
    vector<set<int>> ss_in(K);
    vector<set<int>> ss_out(K);
    double delta = sqrt(-log(1-1.0/N));
    cout << "delta: " << delta << endl;
    for (int i=0; i<N; ++i) {
        if (overlap) {
            for (int k=0; k<K; ++k) {
                if (F(i,k) > delta) {
                    ss[k].insert(i);
                    ss_out[k].insert(i);
                }
                if (H(i,k) > delta) {
                    ss[k].insert(i);
                    ss_in[k].insert(i);
                }
            }
        } else {
            ss_out[index_max(F.row(i))].insert(i);
            ss_in[index_max(H.row(i))].insert(i);
            ss[index_max(F.row(i)+H.row(i))].insert(i);
        }
    }

    ofstream ofs("community_out.txt");
    for (int k=0; k<K; ++k) {
        for (auto itr = ss_out[k].begin(); itr != ss_out[k].end(); ++itr)
            ofs << id2str[*itr] << " ";
        ofs << endl;
    }
    ofs.close();
    ofs.open("community_in.txt");
    for (int k=0; k<K; ++k) {
        for (auto itr = ss_in[k].begin(); itr != ss_in[k].end(); ++itr)
            ofs << id2str[*itr] << " ";
        ofs << endl;
    }
    ofs.close();
    return ss;
}

double CODA::likelihood() const
{
    double lh = .0;
    for (int u=0; u<N; ++u)
        lh += Flikelihood(F.row(u), u);
    return lh;
}
