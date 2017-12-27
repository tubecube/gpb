#include "metrics.hpp"


vector<set<string>> Metrics::file_to_set2(const string& filename)
{
	vector<set<string>> ground;
	ifstream stream(filename);
	if (stream.is_open())
	{
		for (string line; getline(stream, line); )
		{
			int comm;
			string node;
			istringstream fields(line);
			fields >> node >> comm;
			if (ground.size() < comm+1)
				ground.resize(comm+1);
			ground[comm].insert(node);
		}
	} else ERROR("File to set: fail opening file %s\n", filename.c_str());
	return ground;
}


vector<set<string>> Metrics::file_to_set1(const string& filename)
{
	vector<set<string>> ground;
	ifstream stream(filename);
	if (stream.is_open())
	{
		for (string line; getline(stream, line); )
		{
			set<string> s;
			istringstream fields(line);
			for (string field; fields >> field; )
			{
				if (field.length() == 0)
					WARN("File to set: two consecutive tabs, or tab at the start of a line. Ignoring empty fields like this.\n");
				else
					s.insert(field);
			}
			if (s.size() == 0)
				WARN("File to set: ignoring empty sets.\n");
			else
				ground.push_back(s);
		}
	} else ERROR("File to set: fail opening file! %s\n", filename.c_str());

    return ground;
}


void Metrics::set_to_file1(const string& filename, const vector<set<string>>& comm)
{
	ofstream stream(filename);
	for (const set<string>& c : comm)
	{
		for (const string& node : c)
			stream << node << '\t';
		stream << '\n';
	}
	stream.close();
}

void Metrics::set_to_file2(const string& filename, const vector<set<string>>& comm)
{
	ofstream stream(filename);
	for (int i = 0; i < comm.size(); i++)
	{
		for (const string& node : comm[i])
			stream << node << '\t' << i;
		stream << '\n';
	}
	stream.close();
}

void Metrics::keepOnlyOverlap(vector<set<string>>& g1, vector<set<string>>& g2)
{
    std::set<string> set1;
    std::set<string> set2;
    for (int i=g2.size()-1; i>=0; i--)
        set2.insert(g2[i].begin(), g2[i].end());

    for (int i=g1.size()-1; i>=0; i--)
    {
        auto it = g1[i].begin();
        auto end = g1[i].end();
        while (it != end) {
            if (set2.find(*it) == set2.end())
                it = g1[i].erase(it);
            else
                ++it;
        }
    }

    for (int i=g1.size()-1; i>=0; i--)
        set1.insert(g1[i].begin(), g1[i].end());

    for (int i=g2.size()-1; i>=0; i--)
    {
        auto it = g2[i].begin();
        auto end = g2[i].end();
        while (it != end) {
            if (set1.find(*it) == set1.end())
                it = g2[i].erase(it);
            else
                ++it;
        }
    }
}

double Metrics::NNMI(const vector<set<string>>& g1, const vector<set<string>>& g2)
{
    if (is_overlap(g1) || is_overlap(g2))
        WARN("NNMI: input sets may overlap!\n");

    size_t n1 = g1.size();
    size_t n2 = g2.size();
    size_t confusion[n1][n2];

    for (int i=0; i<n1; ++i)
        for (int j=0; j<n2; ++j)
            confusion[i][j] = overlap_size(g1[i], g2[j]);

    size_t N1[n1], N2[n2];
    for (int i=0; i<n1; ++i)
        N1[i] = g1[i].size();
    for (int i=0; i<n2; ++i)
        N2[i] = g2[i].size();

    set<string> ss;
    for (int i=g1.size()-1; i>=0; i--)
        ss.insert(g1[i].begin(), g1[i].end());
    for (int i=g2.size()-1; i>=0; i--)
        ss.insert(g2[i].begin(), g2[i].end());
    size_t N = ss.size();

    double MI = .0;
    double normalizer = .0;

    for (int i=0; i<n1; ++i)
    {
        for (int j=0; j<n2; ++j)
        {
            if (confusion[i][j] == 0)
                continue;
            MI += confusion[i][j] * log((double)confusion[i][j] * N / N1[i] / N2[j]);
        }
    }
    for (int i=0; i<n1; ++i)
    {
        if (N1[i] == 0)
            continue;
        normalizer -= N1[i] * log((double)N1[i] / N);
    }
    for (int j=0; j<n2; ++j)
    {
        if (N2[j] == 0)
            continue;
        normalizer -= N2[j] * log((double)N2[j] / N);
    }
    normalizer /= 2;

    double NMI = MI/normalizer;
    INFO("Metric: NMI: %lf\n",NMI);
    return NMI;
}

vector<double> Metrics::ONMI(const vector<set<string>>& g1, const vector<set<string>>& g2)
{
    if (!is_overlap(g1) && !is_overlap(g2))
        WARN("ONMI: input sets are both disjoint!\n");

    set<string> ss;
    for (int i=g1.size()-1; i>=0; i--)
        ss.insert(g1[i].begin(), g1[i].end());
    for (int i=g2.size()-1; i>=0; i--)
        ss.insert(g2[i].begin(), g2[i].end());
    size_t N = ss.size();

    double lfkIXY, lfkIYX;
    double HX, HY, HXgivenY, HYgivenX;

    onmi_one_turn(g1, g2, N, HX, HXgivenY, lfkIXY);
    onmi_one_turn(g2, g1, N, HY, HYgivenX, lfkIYX);

    double NMILFK = 0.5 * ( lfkIXY + lfkIYX );
    double NMIMAX = 0.5 * ( (HX-HXgivenY)/max(HX, HY) + (HY-HYgivenX)/max(HX, HY) );
    double NMISUM = 0.5 * ( (HX-HXgivenY)/((HX+HY)/2) + (HY-HYgivenX)/((HX+HY)/2) );

    INFO("Metric: NMI_lfk: %lf\n", NMILFK);
    INFO("Metric: NMI_max: %lf\n", NMIMAX);
    INFO("Metric: NMI_sum: %lf\n", NMISUM);

    return vector<double>({NMILFK, NMIMAX, NMISUM});
}

void Metrics::onmi_one_turn(const vector<set<string>>& g1, const vector<set<string>>& g2, size_t N,
                        double& HX, double& HXgivenY, double& lfkIXY)
{
    HX = .0; HXgivenY = .0; lfkIXY = .0;

    for (int i=g1.size()-1; i>=0; --i)
    {
        size_t size_x = g1[i].size();
        double Hx = H( (double)size_x/N ) + H( (double)(N-size_x)/N );
        double HxgivenY = Hx;
        for (int j=g2.size()-1; j>=0; --j)
        {
            size_t size_y = g2[j].size();
            size_t overlap = overlap_size(g1[i], g2[j]);
            size_t xminusy = size_x - overlap;
            size_t yminusx = size_y - overlap;
            size_t nxny = N - xminusy - yminusx - overlap;
            double t1 = H( (double)xminusy/N ) + H( (double)yminusx/N );
            double t2 = H( (double)overlap/N ) + H( (double)nxny/N );
            if (t1 > t2)
                continue;
            double Hy = H( (double)size_y/N ) + H( (double)(N-size_y)/N );
            HxgivenY = min(HxgivenY, t1+t2-Hy);
        }
        HX += Hx;
        HXgivenY += HxgivenY;

        lfkIXY += (Hx == .0 ? 1.0 : (Hx-HxgivenY)/Hx);
    }

    lfkIXY /= g1.size();
}

size_t Metrics::overlap_size(const set<string>& setx, const set<string>& sety)
{
    size_t size_x = setx.size();
    size_t size_y = sety.size();
	size_t min_size = min(size_x, size_y);
	if (min_size == 0)	return 0;
    vector<string> tmp(min_size);
    typename vector<string>::iterator it;
    it = set_intersection(setx.begin(), setx.end(), sety.begin(), sety.end(), tmp.begin());
    size_t overlap = it - tmp.begin();
    return overlap;
}

bool Metrics::is_overlap(const vector<set<string>>& g)
{
    set<string> ss;
    size_t n = g.size();
    size_t N = 0;
    for (int i=0; i<n; ++i)
    {
        ss.insert(g[i].begin(), g[i].end());
        N += g[i].size();
        if (ss.size() != N)
            return true;
    }
    return false;
}

inline double Metrics::H(double p)
{
    return p == .0 ? .0 : -p*log2(p);
}

double Metrics::F1(const vector<set<string>>& g1, const vector<set<string>>& g2)
{
	double score = 0.0;
	size_t g1_total = 0, g2_total = 0;
	for (int i=g1.size()-1; i>=0; i--)
		g1_total += g1[i].size();
	for (int i=g2.size()-1; i>=0; i--)
		g2_total += g2[i].size();

	vector<double> store1(g1.size(), 0.0);
	vector<double> store2(g2.size(), 0.0);

	for (int i=g1.size()-1; i>=0; i--)
	{
		for (int j=g2.size()-1; j>=0; j--)
		{
			double F1 = 0.0;
			size_t overlap = overlap_size(g1[i], g2[j]);
			if (overlap != 0)
			{
				double precision = (double) overlap / g2[j].size();
				double recall = (double) overlap / g1[i].size();
				F1 = 2*precision*recall / (precision+recall);
			}
			store1[i] = max(store1[i], F1);
			store2[j] = max(store2[j], F1);
		}
	}
	for (int i=g1.size()-1; i>=0; i--)
		score += store1[i]*((double)g1[i].size()/g1_total);
	for (int i=g2.size()-1; i>=0; i--)
		score += store2[i]*((double)g2[i].size()/g2_total);

	score /= 2;
	INFO("Metric: F1: %lf\n", score);
	return score;
}

double Metrics::F1_one_turn(const vector< set<string> > &g1, const vector< set<string> >& g2)
{
    double F1_avg = 0.0;
	size_t total_size = 0;
	for (int i=g1.size()-1; i>=0; i--)
		total_size += g1[i].size();
    for (int i=g1.size()-1; i>=0; i--)
    {
        double F1 = 0.0;
        for (int j=g2.size()-1; j>=0; j--)
        {
            size_t overlap = overlap_size(g1[i], g2[j]);
			if (overlap != 0)
			{
				double precision = (double) overlap / g2[j].size();
				double recall = (double) overlap / g1[i].size();
				F1 = max(F1, 2*precision*recall/(precision+recall));
			}
        }
        F1_avg += F1*((double)g1[i].size()/total_size);
    }
    return F1_avg;
}
