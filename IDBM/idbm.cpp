#include "IDBM.h"
#include "NMI.h"

pair<int,int>
IDBM::sampling(pair<int,int> nodepr, pair<int,int> oldpr)
{
    int s = nodepr.first;
    int t = nodepr.second;
    int gold = oldpr.first;
    int hold = oldpr.second;
    N[s][gold]--;
    P[t][hold]--;
    Nt[gold]--;
    Pt[hold]--;
    R[gold][hold]--;
    
    double prob[topicnum][topicnum];
    double prev = .0;
    for (int g=0; g<topicnum; ++g)
    {
        for (int h=0; h<topicnum; ++h)
        {
            prev += (R[g][h]+(alpha/topicnum)) * (N[s][g]+P[s][g]+beta) * (N[t][h]+P[t][h]+beta) / (Nt[g]+Pt[g]+nodenum*beta) / (Nt[h]+Pt[h]+nodenum*beta - (g==h?1:0));
            prob[g][h] = prev;
        }
    }
    generator.seed(rd());
    uniform_real_distribution<> dist(0.0, prev);
    double rnd = dist(generator);
    for (int g=0; g<topicnum; ++g)
        for (int h=0; h<topicnum; ++h)
            if (prob[g][h] >= rnd)
            {
                N[s][g]++;
                P[t][h]++;
                Nt[g]++;
                Pt[h]++;
                R[g][h]++;
                return make_pair(g,h);
            }
    exit(1);
}

pair<int,int>
IDBM::sampling2(pair<int,int> nodepr, pair<int,int> oldpr)
{
    int s = nodepr.first;
    int t = nodepr.second;
    int gold = oldpr.first;
    int hold = oldpr.second;
    N[s][gold]--;
    Nt[gold]--;
    R[gold][hold]--;
    double prob[topicnum];
    double prev = .0;
    for (int g=0; g<topicnum; ++g)
    {
        prev += (R[g][hold]+(alpha/topicnum)) * (N[s][g]+P[s][g]+beta) / (Nt[g]+Pt[g]+nodenum*beta) / (Pt[hold]+alpha);
        prob[g] = prev;
    }
    generator.seed(rd());
    uniform_real_distribution<> dist1(0.0, prev);
    double rnd = dist1(generator);
    int gnew;
    for (gnew=0; gnew<topicnum; ++gnew)
        if (prob[gnew] >= rnd)
            break;

    N[s][gnew]++;
    Nt[gnew]++;
    P[t][hold]--;
    Pt[hold]--;

    prev = .0;
    for (int h=0; h<topicnum; ++h)
    {
        prev += (R[gnew][h]+(alpha/topicnum)) * (N[t][h]+P[t][h]+beta) / (Nt[h]+Pt[h]+nodenum*beta) / (Nt[gnew]+alpha);
        prob[h] = prev;
    }

    uniform_real_distribution<> dist2(0.0, prev);
    rnd = dist2(generator);
    int hnew;
    for (hnew=0; hnew<topicnum; ++hnew)
        if (prob[hnew] >= rnd)
            break;

    P[t][hnew]++;
    Pt[hnew]++;
    R[gnew][hnew]++;

    return make_pair(gnew,hnew);
}


void
IDBM::gibbs(int burnin)
{
    init();
    cout << "start Gibbs sampling" << endl;
    for (iter=1; iter<=burnin; iter++)
    {
        int ni = 0;
        for (int i=0; i<nodenum; ++i)
            for (auto ptr=network[i].cbegin(); ptr!=network[i].cend(); ++ptr)
            {
                if (ni%10000 == 0 && ni!=0)
                    cout << setfill('0') << setw(4) << iter << " " << ni << " edges" << flush << '\r';
                Z[ni++] = sampling(make_pair(i,*ptr), Z[ni]);
            }
        cout << setfill('0') << setw(4) << iter << " " << likelihood() << flush << '\r';
    }
}

void
IDBM::save_PSI(string file)
{
    ofstream ofs(file);
    for (int g=0; g<topicnum; ++g)
    {
        double sum = 0;
        for (int h=0; h<topicnum; ++h)
            sum += R[g][h];
        ofs << sum << ":";
        for (int h=0; h<topicnum; ++h)
            ofs << setprecision(4) << fixed << R[g][h]/sum << " ";
        ofs << endl;
    }
    ofs.close();
}
void
IDBM::init()
{
    if (this->initialized == true)
        return;
    // resize and resize
    N.resize(nodenum);
    for (int i=0; i<nodenum; ++i)
        N[i].resize(topicnum);
    P.resize(nodenum);
    for (int i=0; i<nodenum; ++i)
        P[i].resize(topicnum);
    R.resize(topicnum);
    for (int i=0; i<topicnum; ++i)
        R[i].resize(topicnum);
    Nt.resize(topicnum);
    Pt.resize(topicnum);
    Z.resize(edgenum);
    
    // init and init
    generator.seed(rd());
    uniform_int_distribution<> distribution(0,topicnum-1);
    int ni = 0;
    for (int i=0; i<nodenum; ++i)
    {
        for (auto ptr=network[i].cbegin(); ptr!=network[i].cend(); ++ptr)
        {
            int k1 = distribution(generator);
            int k2 = distribution(generator);
            N[i][k1]++;
            P[*ptr][k2]++;
            R[k1][k2]++;
            Nt[k1]++;
            Pt[k2]++;
            Z[ni++] = make_pair(k1,k2);
        }
    }

    // loglikelihood 的初始值
    baselh = topicnum*(lgamma(nodenum*beta)-nodenum*lgamma(beta));

    // 初始化完成
    this->initialized = true;
}

vector<set<string>>
IDBM::link_community() const
{
    vector<set<string>> ss(topicnum, set<string>());
    for (int i=0; i<nodenum; ++i)
    {
        int maxk = 0;
        int max = 0;
        for (int k=0; k<topicnum; ++k)
        {
            if (max < N[i][k]+P[i][k])
            {
                max = N[i][k]+P[i][k];
                maxk = k;
            }       
        }
        ss[maxk].insert(id2str[i]);
    }
    return ss;
}

vector<set<string>>
IDBM::community() const
{
    vector<set<string>> ss(topicnum, set<string>());
    for (int i=0; i<nodenum; ++i)
    {
        auto max = max_element(N[i].begin(), N[i].end());
        int k = max-N[i].begin();
        ss[k].insert(id2str[i]);
    }
    return ss;
}
        
double
IDBM::likelihood() const
{
    double lh = baselh;
    for (int k=0; k<topicnum; ++k)
    {
        lh -= lgamma(Nt[k]+Pt[k]+nodenum*beta);

        for (int i=0; i<nodenum; ++i)
            lh += lgamma(N[i][k]+P[i][k]+beta);
    }
    return lh;
}

void usage()
{
        cout << "Usage: idbm filename #clusters #iters\n";
        cout << "       -g ground truth filename\n";
        cout << "       -a hyperparameter alpha\n";
        cout << "       -b hyperparameter beta\n";
        cout << "       -u undirected graph\n";

        exit(1);
}

int main(int argc, char* argv[])
{
    char *gn = NULL;
    double alpha = 0.1;
    double beta = 0.1;
    bool directed = true;
    if (argc < 4)
    {
        usage();
    }
    
    char *input = argv[1];
    
    int topicnum = atoi(argv[2]);

    int burnin = atoi(argv[3]);

    for (int i=4; i<argc; i++)
    {
        if (strcmp(argv[i], "-g")==0)
            gn = argv[++i];
        else if (strcmp(argv[i], "-a")==0)
            alpha = atof(argv[++i]);
        else if (strcmp(argv[i], "-b")==0)
            beta = atof(argv[++i]);
        else if (strcmp(argv[i], "-u")==0)
            directed = false;
        else
            usage();
    }

    IDBM idbm(topicnum, alpha, beta);
    idbm.read_from_file(input, directed);
    idbm.gibbs(burnin);

    if (gn)
    {
        vector<set<string>> g1 = fileToSet(gn);
        vector<set<string>> g2 = idbm.link_community();
        NNMI(g1, g2);
    }
    idbm.save_PSI("psi.txt");

    return 0;
}
