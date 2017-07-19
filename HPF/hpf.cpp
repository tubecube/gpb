#include "hpf.hpp"

void
HPF::run(int max_iters)
{
    init();
    double lh;
    int iter = 1;
    time_t start_time = time(NULL);
    for (; max_iters == -1 || iter <= max_iters; ++iter) {

        gnshp.fill(_a);
        gnrte = repmat(kshp/krte, K, 1);
        lnshp.fill(_c);
        lnrte = repmat(tshp/trte, K, 1);

        // compute E[ln(x)]
        compute_exp_ln();

        // estimate phi
        for (int i=0, ei=0; i<N; ++i)
            for (auto itr = network[i].cbegin(), end_itr = network[i].cend(); itr!=end_itr; ++itr, ++ei) {
                if (ei % 100000 == 0 && ei != 0)
                    cout << setfill('0') << setw(3) << iter << ": " << ei << " edges" << flush << '\r';
                int j = *itr;

                // vec phi = compute_phi(i,j);

                if (gconv(i) && lconv(j)) {
                    vector<int> inter;
                    auto git = gactive[i].begin();
                    auto gend = gactive[i].end();
                    auto lit = lactive[j].begin();
                    auto lend = lactive[j].end();
                    while (git != gend && lit != lend) {
                        if(*git < *lit)
                            ++git;
                        else if (*git > *lit)
                            ++lit;
                        else {
                            inter.push_back(*git);
                            ++git;
                            ++lit;
                        }
                    }
                    int size = inter.size();
                    if (size > 0) {
                        vec tmp(size);
                        for (int k=0; k<size; ++k)
                            tmp(k) = exp(Elngamma(inter[k],i) + Elnlambda(inter[k],j));
                        tmp = normalise(tmp,1);
                        for (int k=0; k<size; ++k) {
                            gnshp(inter[k], i) += tmp(k);
                            lnshp(inter[k], j) += tmp(k);
                        }
                    }
                } else if (gconv(i) && !lconv(j)) {
                    int size = gactive[i].size();
                    vec tmp(size);
                    auto it = gactive[i].begin();
                    auto end = gactive[i].end();
                    int k = 0;
                    while (it != end) {
                        tmp(k) = exp(Elngamma(*it,i) + Elnlambda(*it,j));
                        ++it;
                        ++k;
                    }
                    tmp = normalise(tmp,1);
                    k = 0;
                    it = gactive[i].begin();
                    while (it != end) {
                        gnshp(*it, i) += tmp(k);
                        lnshp(*it, j) += tmp(k);
                        ++it;
                        ++k;
                    }
                } else if (!gconv(i) && lconv(j)) {
                    int size = lactive[j].size();
                    vec tmp(size);
                    auto it = lactive[j].begin();
                    auto end = lactive[j].end();
                    int k = 0;
                    while (it != end) {
                        tmp(k) = exp(Elngamma(*it,i) + Elnlambda(*it,j));
                        ++it;
                        ++k;
                    }
                    tmp = normalise(tmp,1);
                    k = 0;
                    it = lactive[j].begin();
                    while (it != end) {
                        gnshp(*it, i) += tmp(k);
                        lnshp(*it, j) += tmp(k);
                        ++it;
                        ++k;
                    }
                } else {
                    vec logexp = Elngamma.col(i) + Elnlambda.col(j);
                    double logsum;
                    for (int k=0; k<K; ++k) {
                        if (k == 0)
                            logsum = logexp(k);
                        else if (logsum > logexp(k))
                            logsum = logsum + log(1 + exp(logexp(k) - logsum));
                        else
                            logsum = logexp(k) + log(1 + exp(logsum - logexp(k)));
                    }
                    logexp = exp(logexp - logsum);
                    gnshp.col(i) += logexp;
                    lnshp.col(j) += logexp;
                }
            }

        swap(gnshp, gshp);

        swap(lnshp, lshp);

        // update gamma,lambda
        
        vec lsum = sum(lshp/lrte,1);
        gnrte += repmat(lsum, 1, N);

        swap(gnrte, grte);
        
        vec gsum = sum(gshp/grte,1);
        lnrte += repmat(gsum, 1, N);

        swap(lnrte, lrte);

        // update kappa,tau

        krte = __b + sum(gshp/grte,0);
        trte = __d + sum(lshp/lrte,0);

        if (iter > 10) check_converge();

        double pre = lh;
        lh = heldout.size() ? validation_likelihood() : likelihood();
        cout << setfill('0') << setw(3) << iter << ": " << lh << '(' << (unsigned long)(time(NULL) - start_time) << "s)" << flush << '\r';
        if (iter > 10 && (lh < pre || abs((lh-pre) / lh) < 1e-5))
            break;
    }
    cout << setfill('0') << setw(3) << iter << ": " << lh << '(' << (unsigned long)(time(NULL) - start_time) << "s)" << flush << '\n';
}

void HPF::check_converge() {
    int converged_count = 0;
    for (int i=0; i<N; ++i)
    {
        if (gconv(i)) {
            auto it = gactive[i].begin();
            auto end = gactive[i].end();
            while (it != end) {
                if (gshp(*it,i) <= _a + network[i].size() * 1e-2)
                    it = gactive[i].erase(it);
                else
                    ++it;
            }
            converged_count++;
        } else {
            int active_count = 0;
            for (int k=0; k<K; ++k)
            {
                if (gshp(k,i) > _a + network[i].size() * 1e-2)
                {
                    gactive[i].push_back(k);
                    active_count++;
                }
            }
            if (active_count > K/10)
                gactive[i].clear();
            else {
                converged_count++;
                gconv(i) = 1;
            }
        }
    }

    // cout << "gamma: " << converged_count << endl;

    converged_count = 0;
    for (int i=0; i<N; ++i) {
        if (lconv(i)) {
            auto it = lactive[i].begin();
            auto end = lactive[i].end();
            while (it != end) {
                if (lshp(*it,i) <= _c + network2[i].size() * 1e-2)
                    it = lactive[i].erase(it);
                else
                    ++it;
            }
            converged_count++;
        } else {
            int active_count = 0;
            for (int k=0; k<K; ++k)
            {
                if (lshp(k,i) > _c + network2[i].size() * 1e-2)
                {
                   lactive[i].push_back(k);
                   active_count++;
                }
            }
            if (active_count > K/10)
                lactive[i].clear();
            else {
                lconv[i] = 1;
                converged_count++;
            }
        }
    }
    // cout << "lambda: " << converged_count << endl;
}

void HPF::compute_exp_ln() {
    int ifault;
    for (int i=0; i<N; ++i) {
        if (gconv(i)) {
            for (auto it = gactive[i].begin(); it != gactive[i].end(); ++it)
                Elngamma(*it,i) = digamma(gshp(*it,i),&ifault) - log(grte(*it,i));
        } else {
            for (int k=0; k<K; ++k)
                Elngamma(k,i) = digamma(gshp(k,i),&ifault) - log(grte(k,i));
        }

        if (lconv(i)) {
            for (auto it = lactive[i].begin(); it != lactive[i].end(); ++it)
                Elnlambda(*it,i) = digamma(lshp(*it,i),&ifault) - log(lrte(*it,i));
        } else {
            for (int k=0; k<K; ++k)
                Elnlambda(k,i) = digamma(lshp(k,i),&ifault) - log(lrte(k,i));
        }
    }
}

void HPF::load_model(string dirname) {
    if (dirname.back() != '/')
        dirname += '/';
    gshp.load(dirname+"gshp.dat");
    grte.load(dirname+"grte.dat");
    lshp.load(dirname+"lshp.dat");
    lrte.load(dirname+"lrte.dat");
    tshp.load(dirname+"tshp.dat");
    trte.load(dirname+"trte.dat");
    kshp.load(dirname+"kshp.dat");
    krte.load(dirname+"krte.dat");
}

void HPF::save_model(string dirname) const {
    if (dirname.back() != '/') 
        dirname = dirname + '/';
    mat tmp;
    tmp = gshp.t();
    tmp.save(dirname+"gshp.dat", arma_ascii);
    tmp = grte.t();
    tmp.save(dirname+"grte.dat", arma_ascii);
    tmp = lshp.t();
    tmp.save(dirname+"lshp.dat", arma_ascii);
    tmp = lrte.t();
    tmp.save(dirname+"lrte.dat", arma_ascii);
    tmp = tshp.t();
    tmp.save(dirname+"tshp.dat", arma_ascii);
    tmp = trte.t();
    tmp.save(dirname+"trte.dat", arma_ascii);
    tmp = kshp.t();
    tmp.save(dirname+"kshp.dat", arma_ascii);
    tmp = krte.t();
    tmp.save(dirname+"krte.dat", arma_ascii);
}

vec HPF::compute_phi(int i, int j) {
    vec logexp = Elngamma.col(i) + Elnlambda.col(j);
    if (gconv(i) && lconv(j)) {
        vector<int> inter;
        auto git = gactive[i].begin();
        auto gend = gactive[i].end();
        auto lit = lactive[j].begin();
        auto lend = lactive[j].end();
        while (git != gend && lit != lend) {
            if(*git < *lit)
                ++git;
            else if (*git > *lit)
                ++lit;
            else {
                inter.push_back(*git);
                ++git;
                ++lit;
            }
        }
        int size = inter.size();
        if (size == 0)
            logexp.zeros();
        else {
            vec tmp(size);
            for (int k=0; k<size; ++k)
                tmp(k) = exp(Elngamma(inter[k],i) + Elnlambda(inter[k],j));
            tmp = normalise(tmp,1);
            logexp.zeros();
            for (int k=0; k<size; ++k)
                logexp(inter[k]) = tmp(k); 
        }
    } else if (gconv(i) && !lconv(j)) {
        int size = gactive[i].size();
        vec tmp(size);
        auto it = gactive[i].begin();
        auto end = gactive[i].end();
        int k = 0;
        while (it != end) {
            tmp(k) = exp(logexp(*it));
            ++it;
            ++k;
        }
        tmp = normalise(tmp,1);
        logexp.zeros();
        k = 0;
        it = gactive[i].begin();
        while (it != end) {
            logexp(*it) = tmp(k);
            ++it;
            ++k;
        }
    } else if (!gconv(i) && lconv(j)) {
        int size = lactive[j].size();
        vec tmp(size);
        auto it = lactive[j].begin();
        auto end = lactive[j].end();
        int k = 0;
        while (it != end) {
            tmp(k) = exp(logexp(*it));
            ++it;
            ++k;
        }
        tmp = normalise(tmp,1);
        logexp.zeros();
        k = 0;
        it = lactive[j].begin();
        while (it != end) {
            logexp(*it) = tmp(k);
            ++it;
            ++k;
        }
    } else {
        double logsum;
        for (int k=0; k<K; ++k) {
            if (k == 0)
                logsum = logexp(k);
            else if (logsum > logexp(k))
                logsum = logsum + log(1 + exp(logexp(k) - logsum));
            else
                logsum = logexp(k) + log(1 + exp(logsum - logexp(k)));
        }
        logexp = exp(logexp - logsum);
    }
    return logexp; 
}

double
HPF::validation_likelihood() const
{
    double d1 = .0;
    double d2 = .0;
    size_t size1 = 0;
    size_t size2 = 0;
    for (auto itr=heldout.begin(); itr!=heldout.end(); ++itr) {
        int i = itr->first.first;
        int j = itr->first.second;
        
        double mean = dot( (gshp.col(i)/grte.col(i)) , (lshp.col(j)/lrte.col(j)) );

        if (itr->second==true) {
            d1 += log(mean)-mean;
            size1++;
        } else {
            d2 -= mean;
            size2++;
        }
    }
    return (d1+d2)/(size1+size2);
}

double
HPF::likelihood() const
{
    double s = .0;
    mat gexp = gshp/grte;
    mat lexp = lshp/lrte;
    for (int i=0; i<N; ++i)
    {
        for (auto itr = network[i].begin(); itr != network[i].end(); ++itr)
        {
            double mean = dot( gexp.col(i) , lexp.col(*itr) );
            s += log(mean);
        }
    }
    s -= dot( sum(gexp,1) , sum(lexp,1) );
    return s;
}

void
HPF::init()
{
    arma_rng::set_seed_random();

    kshp = ones<rowvec>(N) * (K*_a) + __a;
    krte = randu<rowvec>(N) * 0.01 + __b;

    tshp = ones<rowvec>(N) * (K*_c) + __c;
    trte = randu<rowvec>(N) * 0.01 + __d;

    lshp = randu<mat>(K,N) * 0.01 + _c;
    lrte = randu<mat>(K,N) * 0.01 + 1.0;

    gshp = randu<mat>(K,N) * 0.01 + _a;
    grte = randu<mat>(K,N) * 0.01 + 1.0;
    
    lnshp.zeros(K,N);
    lnrte.zeros(K,N);
    gnshp.zeros(K,N);
    gnrte.zeros(K,N);

    Elngamma.zeros(K,N);
    Elnlambda.zeros(K,N);

    gactive.resize(N);
    lactive.resize(N);

    gconv = zeros<uvec>(N);
    lconv = zeros<uvec>(N);
}

vector<set<int>>
HPF::node_community(bool overlap, int type)
{
    vector<set<int>> ss(K);
    vector<set<int>> ss_in(K);
    vector<set<int>> ss_out(K);
    double delta = sqrt(-log(1-1.0/N));
    if (overlap) {
        for (int i=0; i<N; ++i) {
            if (ground_set.find(i) != ground_set.end()) {
                for (int k=0; k<K; ++k) {
                    if (lshp(k,i)/lrte(k,i) > delta) {
                        ss[k].insert(i);
                        ss_in[k].insert(i);
                    }
                    if (gshp(k,i)/grte(k,i) > delta) {
                        ss[k].insert(i);
                        ss_out[k].insert(i);
                    }
                }
            }
        }
    } else {
        urowvec index = index_max(lshp/lrte + gshp/grte, 0); 
        urowvec index_in = index_max(lshp/lrte, 0);
        urowvec index_out = index_max(gshp/grte, 0);
        for (int i=0; i<N; ++i)
            if (ground_set.find(i) != ground_set.end()) {
                ss[index(i)].insert(i);
                ss_in[index_in(i)].insert(i);
                ss_out[index_out(i)].insert(i);
            }
    }

    if (type == 0)
        return ss;
    else if (type == 1)
        return ss_out;
    else if (type == 2)
        return ss_in;
    else
        exit(1);
}

vector<set<int>>
HPF::link_community(bool overlap, double thresh)
{
    vector<set<int>> ss(K);
    vector<set<int>> ss_in(K);
    vector<set<int>> ss_out(K);
    mat component_in = zeros<mat>(K,N);
    mat component_out = zeros<mat>(K,N);
    mat component = zeros<mat>(K,N);

    for (int i=0; i<N; ++i) {
        for (auto itr=network[i].cbegin(); itr!=network[i].cend(); ++itr) {
            vec phi = compute_phi(i, *itr);
            component_out.col(i) += phi;
            component_in.col(*itr) += phi;
            component.col(i) += phi;
            component.col(*itr) += phi;
        }
    }

    component_out = normalise(component_out, 1, 0);
    component_in = normalise(component_in, 1, 0);
    component = normalise(component, 1, 0);

    // string filename = "component" + to_string(K);
    // ofstream ofs(filename);

    if (overlap) {
        for (int i=0; i<N; ++i)
            for (int k=0; k<K; ++k) {
                if (component_out(k,i) > thresh)
                    ss_out[k].insert(i);
                if (component_in(k,i) > thresh)
                    ss_in[k].insert(i);
                if (component(k,i) > thresh)
                    ss[k].insert(i);
            }
    } else {
        urowvec index = index_max(component, 0);
        urowvec index_in = index_max(component_in, 0);
        urowvec index_out = index_max(component_out, 0);
        for (int i=0; i<N; ++i)
        {
            ss[index(i)].insert(i);
            ss_in[index_in(i)].insert(i);
            ss_out[index_out(i)].insert(i);
        }
    }

    return ss;
}
