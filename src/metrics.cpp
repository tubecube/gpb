template <class T>
vector<set<T>> Metrics<T>::file_to_set2(const string& filename)
{
	vector<set<T>> ground;
    ifstream stream(filename);
	if (stream.is_open())
	{
		for (string line; getline(stream, line); )
		{
			int comm;
			T node;
			istringstream fields(line);
			fields >> node >> comm;
			if (ground.size() < comm+1)
				ground.resize(comm+1);
			ground[comm].insert(node);
		}
	} else ERROR("Metric: fail opening file! [%s]\n", filename.c_str());
	return ground;
}


template <class T>
vector<set<T>> Metrics<T>::file_to_set1(const string& filename)
{
	vector<set<T>> ground;
    ifstream stream(filename);
	if (stream.is_open())
	{
		for (string line; getline(stream, line); )
		{
			set<T> s;
			istringstream fields(line);
			for (T field; fields >> field; )
			{
				if (field.length() == 0)
					WARN("Metric: two consecutive tabs, or tab at the start of a line. Ignoring empty fields like this.\n");
				else
					s.insert(field);
			}
			if (s.size() == 0)
				WARN("Metric: ignoring empty sets.\n");
			else
				ground.push_back(s);
		}
	} else ERROR("Metric: fail opening file! [%s]\n", filename.c_str());

    return ground;
}

template <class T>
void Metrics<T>::set_to_file(const string& filename, const vector<set<T>>& comm)
{
	ofstream stream(filename);
	for (const set<T>& c : comm)
	{
		for (const T& node : c)
			stream << node << '\t';
		stream << '\n';
	}
	stream.close();
}

template <class T>
void Metrics<T>::keepOnlyOverlap(vector<set<T>>& g1, vector<set<T>>& g2)
{
    std::set<T> set1;
    std::set<T> set2;
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

template <class T>
double Metrics<T>::NNMI(const vector< set<T> >& g1, const vector< set<T> >& g2)
{
    if (is_overlap(g1) || is_overlap(g2))
        WARN("Metric: group1 or group2 may overlap!\n");

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

    set<T> ss;
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

template<class T>
double Metrics<T>::ONMI(const vector< set<T> >& g1, const vector< set<T> >& g2)
{
    if (!is_overlap(g1) && !is_overlap(g2))
        WARN("Metric: group1 and group2 are both disjoint!\n");
    
    set<T> ss;
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

    return NMILFK;
}

template <class T>
void Metrics<T>::onmi_one_turn(const vector< set<T> >& g1, const vector< set<T> >& g2, size_t N,
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

template <class T>
size_t Metrics<T>::overlap_size(const set<T>& setx, const set<T>& sety)
{
    size_t size_x = setx.size();
    size_t size_y = sety.size();
    vector<T> tmp( min( size_x, size_y ) );
    typename vector<T>::iterator it;
    it = set_intersection(setx.begin(), setx.end(), sety.begin(), sety.end(), tmp.begin());
    size_t overlap = it - tmp.begin();
    return overlap;
}

template<class T>
bool Metrics<T>::is_overlap(const vector< set<T> >& g)
{
    set<T> ss;
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

template<class T>
inline double Metrics<T>::H(double p)
{
    return p == .0 ? .0 : -p*log2(p);
}

template <class T>
double Metrics<T>::F1(const vector< set<T> >& g1, const vector< set<T> >& g2)
{
    double score = 0.5 * (F1_one_turn(g1, g2) + F1_one_turn(g2, g1));
    INFO("Metric: F1: %lf\n", score);
    return score;
}

template <class T>
double Metrics<T>::F1_one_turn(const vector< set<T> > &g1, const vector< set<T> >& g2)
{
    double F1_avg = 0.0; 
    for (int i=g1.size()-1; i>=0; i--)
    {
        double F1 = 0.0;
        for (int j=g2.size()-1; j>=0; j--)
        {
            size_t overlap = overlap_size(g1[i], g2[j]);
            double precision = (double) overlap / g2[j].size();
            double recall = (double) overlap / g1[i].size();
            F1 = max(F1, 2*precision*recall/(precision+recall));
        }
        F1_avg += F1;
    }
    return F1_avg / g1.size();
}
