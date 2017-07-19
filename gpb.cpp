#include "gpb.hpp"

void GPB::mcmc(int burnin, int sample)
{
    init();
    double lh;
    cout << "burnin: " << burnin << "   sample: " << sample << endl;
    for (iter=1; iter<=burnin; iter++) { 
        Bshp.fill(Bpshp);
        Fshp.fill(Fpshp);
        Brte.fill(Bprte);
        Frte.fill(Fprte);
        for (int i=0; i<N; ++i)
        {
            for (auto it = network[i].cbegin(); it != network[i].cend(); ++it)
            {
                int j = *it;
                mat phi = sample_phi(i, j);
                Fshp.col(i) += sum(phi, 1);
                Fshp.col(j) += sum(phi, 0).t();
                Bshp += phi;
            }
        }
        Frte += (B+B.t())*(repmat(sum(F,1), 1, N) - F);
        F = sample_F(Fshp, Frte);
        EF = F;
        Brte += sum(F,1)*sum(F,1).t() - F*F.t();
        if (iter <= burnin/2)
            Brte = 2 * Brte;
        B = sample_B(Bshp, Brte);
        EB = B;
        lh = heldout.size() ? validation_likelihood() : likelihood();
        cout << setfill('0') << setw(3) << iter << ": " << lh << flush << '\n';
    }

    for (int s=1; s<=sample; s++) {
        Bshp.fill(Bpshp);
        Fshp.fill(Fpshp);
        Brte.fill(Bprte);
        Frte.fill(Fprte);
        for (int i=0; i<N; ++i)
        {
            for (auto it = network[i].cbegin(); it != network[i].cend(); ++it)
            {
                int j = *it;
                mat phi = sample_phi(i, j);
                Fshp.col(i) += sum(phi, 1);
                Fshp.col(j) += sum(phi, 0).t();
                Bshp += phi;
            }
        }
        Frte += (B+B.t())*(repmat(sum(F,1), 1, N) - F);
        F = sample_F(Fshp, Frte);
        Brte += sum(F,1)*sum(F,1).t() - F*F.t();
        B = sample_B(Bshp, Brte);
        EF = EF*((s-1.0)/s) + F/s;
        EB = EB*((s-1.0)/s) + B/s;
        lh = heldout.size() ? validation_likelihood() : likelihood();
        cout << setfill('0') << setw(3) << s << ": " << lh << flush << '\n';
    }
}

mat GPB::sample_phi(int i, int j)
{
    mat phi = zeros<mat>(K,K);
    mat factor = (F.col(i) * F.col(j).t()) % B;
    double sum = accu(factor);
    // cout << i << " " << j << " " << sum << endl;
    poisson_distribution<int> pd(sum);
    int num;
    if (sum < 1) {
        uniform_real_distribution<double> ud(0, 1);
        do {
            num = pd(generator) + 1;
        } while (ud(generator) >= 1.0/num);
    } else { 
        do {
            num = pd(generator);
        } while (num == 0);
    }
    /*
    uniform_real_distribution<double> urd(0.0, sum);
    while (num > 0) {
        double s = urd(generator);
        double acc = 0;
        for (int k1=0; k1<K; ++k1) {
            for (int k2=0; k2<K; ++k2) {
                acc += factor(k1,k2);
                if (acc > s) {
                    phi(k1,k2) += 1;
                    goto end;
                }
            }
        }
end: num--;
    }
    */
    phi = (factor/sum)*num;
    return phi;
}

mat GPB::sample_F(mat Fshp, mat Frte)
{
    mat F = zeros<mat>(K,N);
    generator.seed(rd());
    for (int n=0; n<N; ++n) {
        for (int k=0; k<K; ++k) {
            gamma_distribution<double> gd(Fshp(k,n), 1.0/Frte(k,n));
            F(k,n) = gd(generator);
        }
    }
    return F;
}

mat GPB::sample_B(mat Bshp, mat Brte)
{
    // cout << Bshp;
    // cout << Brte;
    mat B = zeros<mat>(K,K);
    generator.seed(rd());
    for (int k1=0; k1<K; ++k1) {
        for (int k2=0; k2<K; ++k2) {
            gamma_distribution<double> gd(Bshp(k1,k2), 1.0/Brte(k1,k2));
            B(k1,k2) = gd(generator);
        }
    }
    // cout << B;
    return B;
}

void GPB::run(int max_iters)
{
    init();
    double lh;
    time_t start_time = time(NULL);
    time_t last_time = start_time;
    mat oldEF;
    double delta;
    double li;
    for (iter=1; max_iters == -1 || iter <= max_iters; ++iter) {
        do {
            Fshp.fill(Fpshp);
            Bshp.fill(Bpshp);
            for (int i=0, ei=1; i<N; ++i) {
                for (auto it = network[i].cbegin(), end = network[i].cend(); it != end; ++it, ++ei) {
                    if (ei % 100000 == 0)
                        cout << setfill('0') << setw(3) << iter << ": " << ei << " edges" << flush << '\r';
                    int j = *it;
                    
                    mat phi = estimate_phi(i,j);
                    Bshp += phi;
                    // cout << phi;
                    Fshp.col(i) += sum(phi,1);
                    Fshp.col(j) += sum(phi,0).t();
                }
            }

            Frte.fill(Fprte);
            Frte += (EB+EB.t())*(repmat(sum(EF,1), 1, N) - EF);
            oldEF = move(EF);
            EF = Fshp/Frte;
            li = heldout.size() ? validation_likelihood() : likelihood();
            // cout << setfill('0') << setw(3) << iter << ": " << li << '(' << (time(NULL) - last_time) << "s)" << endl;
            delta = norm(oldEF-EF, 1);
            // cout << delta << endl;
        } while (delta >= 1);

        Brte.fill(Bprte);
        Brte += sum(EF,1)*sum(EF,1).t() - EF*EF.t();
        EB = Bshp/Brte;
        // cout << EB;

        compute_exp_ln();

        double pre = lh;
        lh = heldout.size() ? validation_likelihood() : likelihood();
        cout << setfill('0') << setw(3) << iter << ": " << lh << '(' << (time(NULL) - last_time) << "s)" << endl;
        last_time = time(NULL);
        if (iter > 10 && (lh < pre || abs((lh-pre) / lh) < 1e-5))
            break;
    }
    cout << "total seconds: " << last_time - start_time << endl;
}

mat GPB::estimate_phi(int i, int j) {
    mat phi(K,K);
    double logsum;
    bool first = true;
    for (int k1=0; k1<K; ++k1) {
        for (int k2=0; k2<K; ++k2) {
            phi(k1,k2) = ElnF(k1,i) + ElnF(k2,j) + ElnB(k1,k2);
            if (first) {
                logsum = phi(k1,k2);
                first = false;
            } else if (phi(k1,k2) < logsum) {
                logsum = logsum + log(1 + exp(phi(k1,k2) - logsum));
            } else {
                logsum = phi(k1,k2) + log(1 + exp(logsum - phi(k1,k2)));
            }
        }
    }
    return exp(phi-logsum);
}
/*
sp_mat GPB::estimate_phi(int i, int j) {
    sp_mat phi(K,K);
    double logsum;
    sp_mat::iterator end = B.end();
    bool first = true;
    for (sp_mat::iterator it = B.begin(); it!=end; ++it) {
        int k1 = it.row();
        int k2 = it.col();

        phi(k1,k2) = ElnF(k1,i) + ElnF(k2,j) + ElnB(k1,k2);
        if (first) {
            logsum = phi(k1,k2);
            first = false;
        }
        else if (phi(k1,k2) < logsum)
            logsum = logsum + log(1 + exp(phi(k1,k2) - logsum));
        else
            logsum = phi(k1,k2) + log(1 + exp(logsum - phi(k1,k2)));
    }
    for (sp_mat::iterator it = B.begin(); it!=end; ++it) {
        int k1 = it.row();
        int k2 = it.col();
        phi(k1,k2) = exp(phi(k1,k2)-logsum);
    }
    return phi;
}
*/
void GPB::save_model(string dirname) const {
    if (dirname.back() != '/') 
        dirname = dirname + '/';
    EF.save(dirname + "F.dat", arma_ascii);
    EB.save(dirname + "B.dat", arma_ascii);
}

void GPB::init() {
    generator.seed(rd());

    arma_rng::set_seed_random();

    Fshp = ones<mat>(K,N) * Fpshp;
    Frte = ones<mat>(K,N) * Fprte;
    EF = Fshp/Frte;
    // F = sample_F(Fshp, Frte);
    F = EF;
    ElnF.set_size(K,N);

    // Bshp = eye<mat>(K,K) * 0.8 + randu<mat>(K,K) * 0.01 + Bpshp;
    Bshp = ones<mat>(K,K) * Bpshp;
    Brte = ones<mat>(K,K) * Bprte;
    EB = Bshp/Brte;
    // B = sample_B(Bshp, Brte);
    B = EB;
    // B = eye<mat>(K,K) * 0.01;
    ElnB.set_size(K,K);

    compute_exp_ln();

    //--------------------//
    umat location(2,M);
    for (int i=0, ei=0; i<N; ++i) {
        for (auto it=network[i].begin(); it!=network[i].end(); ++it, ++ei) {
            location(0,ei) = i;
            location(1,ei) = *it;
        }
    }
    ivec values = ones<ivec>(M);
    sp_imat X(location,values,N,N);
    this->Network = X;
    //-------------------//
}

double GPB::validation_likelihood() const {
    double d1 = .0;
    double d2 = .0;
    size_t size1 = 0;
    size_t size2 = 0;
    for (auto itr = heldout.cbegin(); itr!=heldout.cend(); ++itr) {
        int i = itr->first.first;
        int j = itr->first.second;

        double mean = accu((EF.col(i)).t() * EB * EF.col(j));

        if (itr->second) {
            d1 += log(mean) - mean;
            size1++;
        } else {
            d2 -= mean;
            size2++;
        }
    }
    return (d1+d2)/(size1+size2);
}

double GPB::likelihood() const {
    double s = .0;
    for (int i=0; i<N; ++i)
    {
        for (auto itr = network[i].begin(); itr != network[i].end(); ++itr)
        {
            int j = *itr;
            double mean = accu((EF.col(i)).t() * EB * EF.col(j));
            s += log(mean);
        }
    }
    s -= accu(sum(EF,1).t() * EB * sum(EF,1));
    return s;
}

void GPB::compute_exp_ln() {
    int ifault;
    for (int i=0; i<N; ++i)
        for (int k=0; k<K; ++k)
            ElnF(k,i) = digamma(Fshp(k,i),&ifault) - log(Frte(k,i));

    for (int k1=0; k1<K; ++k1)
        for (int k2=0; k2<K; ++k2)
            ElnB(k1,k2) = digamma(Bshp(k1,k2),&ifault) - log(Brte(k1,k2));
}

vector<set<int>> GPB::link_community(bool overlap, double thresh)
{
    vector<set<int>> ss(K);
    mat component = zeros<mat>(K,N);

    for (int i=0; i<N; ++i) {
        for (auto it=network[i].cbegin(); it!=network[i].cend(); ++it) {
            mat phi = sample_phi(i, *it);
            // cout << phi;
            component.col(i) += sum(phi,1);
            component.col(*it) += sum(phi,0).t();
        }
    }

    if (overlap) {
        for (int i=0; i<N; ++i) {
            for (int k=0; k<K; ++k)
                if (component(k,i) > thresh)
                    ss[k].insert(i);
        }
    } else {
        urowvec index = index_max(component, 0);
        for (int i=0; i<N; ++i)
            ss[index(i)].insert(i);
    }
    return ss;
}

vector<vector<int>> GPB::node_community(bool overlap, double thresh)
{
    if (thresh == 0)
        thresh = sqrt(-log(1 - 1.0/N));
    vector<vector<int>> ss(K);
    // mat component = normalise(EF, 1, 0);
    // cout << component;
    if (overlap) {
        
        for (int i=0; i<N; ++i)
            for (int k=0; k<K; ++k)
                if (EF(k,i) > thresh)
                    ss[k].push_back(i);
        /*
        for (int k=0; k<K; ++k) {
            vec cur = EF.row(k).t();
            uvec indices = sort_index(cur, "descend");
            for (int i=0; i<N; ++i)
                ss[k].push_back(indices(i));
        */
    } else {
        urowvec index = index_max(EF, 0);
        for (int i=0; i<N; ++i)
            ss[index(i)].push_back(i);
    }
    return ss;
}

void GPB::save_community(vector<vector<int>>& comm) {
    ofstream ofs("communities.txt");
    for (int k=0; k<K; k++) {
        for (int idx=0; idx<comm[k].size(); idx++)
            ofs << id2str[comm[k][idx]] << ":" << EF(k,comm[k][idx]) << "\t";
        ofs << endl;
    }
    ofs.close();
}
