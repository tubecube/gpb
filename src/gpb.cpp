#include "gpb.hpp"

int GPB::ECHO_PER_ITERS = 5;

int GPB::SAVE_PER_ITERS = 20;

void GPB::gibbs(int burnin, int Ns)
{
	INFO("GPB: start Gibbs sampling burnin: %d samples: %d\n", burnin, Ns);

	mat EF(size(F), fill::zeros);
	mat EB(size(B), fill::zeros);

	double final_elbo = 0;

	bool directed = graph.is_directed();

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
					save("sample_"+to_string(s), EF, EB);
				}
			}
			if ((iter-1) % ECHO_PER_ITERS == 0)
				INFO("GPB: likelihood after %d iters: %lf\n", iter-1, elbo);
		}

        Bshp.fill(Bpshp);
        Fshp.fill(Fpshp);
        Brte.fill(Bprte);
        Frte.fill(Fprte);

		for (auto it = graph.network.begin(); it != graph.network.end(); ++it)
		{
        	int i = it.row();
			int j = it.col();
			// omit heldouts
			if (graph.check_in_heldouts(i, j, true))
				continue;
            sample_phi(i, j, sample_deep, false);
            Fshp.col(i) += sum(phi, 1);
            Fshp.col(j) += sum(phi, 0).t();
			if (directed)
            	Bshp += phi;
			else
				Bshp += (phi+phi.t()) / 2;
		}
		
		// additional samples for zeros
		if (beta != 0.0)
			graph.push_heldout_with_percentage(beta, 0.0);
		// additional samples in heldouts
		for (const Graph::Heldout& heldout : graph.heldouts)
		{
			for (int i = 0; i < 2; i++)
			{
				for (const pair<int,int>& pr : heldout.pairs[i])
				{
					int source = pr.first;
					int dest = pr.second;
					sample_phi(source,dest,sample_deep,true);
					if (accu(phi) == 0)
					{
						continue;
					}
					Fshp.col(source) += sum(phi, 1);
					Fshp.col(dest) += sum(phi, 0).t();
					if (directed)
						Bshp += phi;
					else
						Bshp += (phi+phi.t()) / 2;
				}
			}
		}
		if (beta != 0.0)
			graph.pop_heldout();


		if (directed)
        	Frte += (B+B.t())*(repmat(sum(F,1), 1, N) - F);
		else
        	Frte += B*(repmat(sum(F,1), 1, N) - F);
        sample_F(Fshp, Frte);

		if (directed)
			Brte += sum(F,1)*sum(F,1).t() - F*F.t();
		else
			Brte += (sum(F,1)*sum(F,1).t() - F*F.t()) / 2;
        sample_B(Bshp, Brte);
    }
	F = EF; B = EB;
	INFO("GPB: final likelihood after %d samples: %lf\n", Ns, final_elbo);

	save("sample_final", EF, EB);
}

void GPB::sample_phi(int i, int j, bool sample_deep, bool accept_zero)
{
	phi.zeros();
    mat factors = (F.col(i) * F.col(j).t()) % B;
    double rate = accu(factors);
	poisson_distribution<int> pd(rate);
	uniform_real_distribution<double> ud(0,1);

    int num;
	if (accept_zero)
	{
		num = pd(generator);
		if (num == 0)
			return;
	}
	else
	{
		if (rate < 1)
		{
			do num = pd(generator) + 1;
			while (ud(generator) * num >= 1.0);
		}
		else
		{ 
			do num = pd(generator);
			while (num == 0);
		}
	}

	factors = factors/rate; // normalize
	if (sample_deep)
	{
		int n_elem = factors.n_elem;
		vector<double> stairs(n_elem);
		stairs[0] = factors(0);
		for (int i = 1; i < n_elem; ++i)
			stairs[i] = stairs[i-1] + factors(i);
		while (num-- > 0)
		{
			vector<double>::iterator iter = std::lower_bound(stairs.begin(), stairs.end(), ud(generator));
			phi[iter-stairs.begin()] += 1;
		}
	}
	else
		phi = factors*num;
}

void GPB::sample_F(const mat& Fshp, const mat& Frte)
{
	int rows = Fshp.n_rows;
	int cols = Fshp.n_cols;
	F.zeros();
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
	B.zeros();
	bool directed = graph.is_directed();
	for (int k1=0; k1<rows; ++k1)
	{
		for (int k2 = directed ? 0 : k1; k2<cols; ++k2)
		{
			gamma_distribution<double> gd(Bshp(k1,k2), 1.0/Brte(k1,k2));
			B(k1,k2) = gd(generator);
		}
	}
	if (!directed)
	{
		B = B + B.t();
		B.diag() /= 2;
	}
}

void GPB::init()
{
    set_seed(generator);
    arma_rng::set_seed_random();
	N = graph.n_nodes();

	phi.set_size(K,K);

	F.set_size(K,N);
    Fshp = randu<mat>(K,N) * Fpshp;
    Frte = randu<mat>(K,N) * Fprte;
    sample_F(Fshp, Frte);

	B.set_size(K,K);
    Bshp = randu<mat>(K,K) % ((alpha-1)*eye<mat>(K,K) + ones<mat>(K,K)) * Bpshp;
    Brte = randu<mat>(K,K) * Bprte;
    sample_B(Bshp, Brte);

	save("sample_init", F, B);
}

vector<pair<float,int>> GPB::link_prediction(const Graph::Heldout& test, bool save_in_file, const string& filename) const
{
	DEBUG("GPB: predict links in test set\n");

	vector<pair<float,int>> pool;

	ofstream stream;
	if (save_in_file)
		stream.open(save_dir+"/"+filename);

	for (int i = 0; i < 2; i++)
		for (const pair<int,int>& pr : test.pairs[i])
		{
			int source = pr.first;
			int dest = pr.second;
			float mean = accu(F.col(source).t() * B * F.col(dest));
			float one_prob = 1 - exp(-mean);
			pool.push_back(make_pair(one_prob,i));
			if (save_in_file)
				stream << graph.get_str(source) << '\t' << graph.get_str(dest)
					<< '\t' << i << '\t' << one_prob << endl;
		}

	if (save_in_file)
		stream.close();

	return pool;
}

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
	F.save(save_dir+"/"+prefix+"F.txt", arma_ascii);
	B.save(save_dir+"/"+prefix+"B.txt", arma_ascii);
}

void GPB::load(const string& prefix)
{
	F.load(save_dir+"/"+prefix+"F.txt");
	B.load(save_dir+"/"+prefix+"B.txt");
}
