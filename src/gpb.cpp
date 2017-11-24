#include "gpb.hpp"

int GPB::ECHO_PER_ITERS = 1;

int GPB::SAVE_PER_ITERS = 10;

void GPB::gibbs(int burnin, int Ns)
{
	INFO("Start Gibbs sampling! [burnin: %d samples: %d]\n", burnin, Ns);

	mat EF(size(F), fill::zeros);
	mat EB(size(B), fill::zeros);

	double final_elbo = 0;

    for (int iter=1; iter<=burnin+Ns; iter++)
	{
		if ((iter-1) > burnin || (iter-1) % ECHO_PER_ITERS == 0)
		{
			double elbo = compute_elbo(F, B);
			if (iter-1 > burnin)
			{
				int s = iter- 1 - burnin;
				final_elbo = final_elbo*((s-1.0)/s) + elbo/s;
        		EF = EF*((s-1.0)/s) + F/s;
        		EB = EB*((s-1.0)/s) + B/s;
				if (s % SAVE_PER_ITERS == 0)
				{
					char tmp[50];
					sprintf(tmp, "%s/%d", dir.c_str(), s);
					save(tmp, EF, EB);
				}
			}

			if ((iter-1) % ECHO_PER_ITERS == 0)
				INFO("Likelihood after %d iters: %lf\n", iter-1, elbo);
		}

        Bshp.fill(Bpshp);
        Fshp.fill(Fpshp);
        Brte.fill(Bprte);
        Frte.fill(Fprte);

		for (auto it = graph.network.begin(); it != graph.network.end(); ++it)
		{
        	int i = it.row();
			int j = it.col();
            mat phi = sample_phi(i, j, sample_deep);
            Fshp.col(i) += sum(phi, 1);
            Fshp.col(j) += sum(phi, 0).t();
            Bshp += phi;
		}

        Frte += (B+B.t())*(repmat(sum(F,1), 1, N) - F);
        sample_F(Fshp, Frte);

        Brte += sum(F,1)*sum(F,1).t() - F*F.t();
        sample_B(Bshp, Brte);
    }

	INFO("Gibbs sampling finished, final likelihood: %lf.\n", final_elbo);

	char tmp[50];
	sprintf(tmp, "%s/final", dir.c_str());
	save(tmp, EF, EB);
}

mat GPB::sample_phi(int i, int j, bool sample_deep)
{
    mat factors = (F.col(i) * F.col(j).t()) % B;
    double rate = accu(factors);
	factors = factors/rate; // normalize
	poisson_distribution<int> pd(rate);
	uniform_real_distribution<double> ud(0,1);
    int num;
    if (rate < 1)
	{
        do
		{
            num = pd(generator) + 1;
        } while (ud(generator) * num >= 1.0);
    }
	else
	{ 
        do
		{
            num = pd(generator);
        } while (num == 0);
    }

	if (sample_deep)
	{
		int n_elem = factors.n_elem;
		vector<double> stairs(n_elem);
		stairs[0] = factors(0);
		for (int i = 1; i < n_elem; ++i)
			stairs[i] = stairs[i-1] + factors(i);
		factors.zeros();
		while (num-- > 0)
		{
			auto iter = std::lower_bound(stairs.begin(), stairs.end(), ud(generator));
			factors[iter-stairs.begin()] += 1;
		}
	}
	else
		factors *= num;

    return factors;
}

void GPB::sample_F(const mat& Fshp, const mat& Frte)
{
	int rows = Fshp.n_rows;
	int cols = Fshp.n_cols;
    for (int n=0; n<cols; ++n)
	{
        for (int k=0; k<rows; ++k)
		{
            gamma_distribution<double> gd(Fshp(k,n), 1.0/Frte(k,n));
            F(k,n) = gd(generator);
        }
    }
}

void GPB::sample_B(const mat& Bshp, const mat& Brte)
{
	int rows = Bshp.n_rows;
	int cols = Bshp.n_cols;
    for (int k1=0; k1<rows; ++k1)
	{
        for (int k2=0; k2<cols; ++k2)
		{
            gamma_distribution<double> gd(Bshp(k1,k2), 1.0/Brte(k1,k2));
            B(k1,k2) = gd(generator);
        }
    }
}

/*
void GPB::vi(int max_iters)
{
    double lh;
    time_t start_time = time(NULL);
    time_t last_time = start_time;
    mat oldEF;
    double delta;
    double li;
    for (int iter=1; max_iters == -1 || iter <= max_iters; ++iter) {
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
*/

/*
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
*/
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

void GPB::init(int alpha)
{
    set_seed(generator);

    arma_rng::set_seed_random();

	N = graph.n_nodes();

	F.set_size(K,N);
    Fshp = randu<mat>(K,N) * Fpshp;
    Frte = randu<mat>(K,N) * Fprte;
    sample_F(Fshp, Frte);
    ElnF.set_size(K,N);

	B.set_size(K,K);
    Bshp = randu<mat>(K,K) % ((alpha-1)*eye<mat>(K,K) + ones<mat>(K,K) * Bpshp);
    Brte = randu<mat>(K,K) % (ones<mat>(K,K) * Bprte);
    sample_B(Bshp, Brte);
    ElnB.set_size(K,K);
    // compute_exp_ln();
}

/*
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
*/

double GPB::compute_elbo(const mat& F, const mat& B) const
{
    double s = .0;
	for (auto it=graph.network.begin(); it!=graph.network.end(); ++it)
	{
		int i = it.row();
        int j = it.col();
        double mean = accu( F.col(i).t() * B * F.col(j) );
		double one_prob = 1 - exp(-mean);
        s += log(one_prob);
		s -= (-mean);
	}
	vec tmp = sum(F,1);
	s += -accu(tmp.t() * B * tmp);
    return s;
}

/*
void GPB::compute_exp_ln() {
    int ifault;
    for (int i=0; i<N; ++i)
        for (int k=0; k<K; ++k)
            ElnF(k,i) = digamma(Fshp(k,i),&ifault) - log(Frte(k,i));

    for (int k1=0; k1<K; ++k1)
        for (int k2=0; k2<K; ++k2)
            ElnB(k1,k2) = digamma(Bshp(k1,k2),&ifault) - log(Brte(k1,k2));
}
*/

mat GPB::link_component()
{
    mat component = zeros<mat>(K,N);

	for (auto it=graph.network.begin(); it!=graph.network.end(); ++it)
	{
		int i = it.row();
		int j = it.col();
        mat phi = F.col(i) * F.col(j).t() % B;
        component.col(i) += sum(phi,1);
        component.col(j) += sum(phi,0).t();
    }
	return component;
}

vector<set<string>> GPB::get_community(const mat& component, bool overlap)
{
	vector<set<string>> comm(K); 
	for (int i = 0; i < N; i++)
	{
		if (overlap)
		{
			for (int k=0; k<K; k++)
				if (component(k, i) > 1)
					comm[k].insert(graph.get_str(i));
		} 
		else
		{
			comm[component.col(i).index_max()].insert(graph.get_str(i));
		}
	}
	return comm;
}

void GPB::save(const string& prefix, const mat& F, const mat& B) const
{
	F.save(prefix+"_F.txt", arma_ascii);
	B.save(prefix+"_B.txt", arma_ascii);
}

void GPB::load(const string& prefix)
{
	F.load(prefix+"_F.txt");
	B.load(prefix+"_B.txt");
}
