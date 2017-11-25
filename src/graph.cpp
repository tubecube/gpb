#include "graph.hpp"

int Graph::read_from_file(const string& filename, bool directed)
{
    ifstream stream(filename);
    if (stream.is_open() == false)
	{
       ERROR("Fail opening file! [%s]\n", filename.c_str());
	   return 2;
	}

    this->directed = directed;

    if (directed)
       INFO("Reading from directed graph! [%s]\n", filename.c_str());
    else
       INFO("Reading from undirected graph! [%s]\n", filename.c_str());

	str2id.clear();
	id2str.clear();

	unsigned m = 0;
    for (string line; getline(stream, line); )
    {
		/*add edge*/
        istringstream edge(line);
        string source, dest;
        edge >> source >> dest;

		if (source.size()==0 || dest.size()==0)
			continue;
		// comment
		if (source[0]=='%' || dest[0]=='%')
			continue;

		int sid = check_id(source);
		int did = check_id(dest);

		N = str2id.size();
		if (edges.size() < N)
			edges.resize(2*N);

		if (!directed && sid > did)
			std::swap(sid, did);
		edges[sid].insert(did);
		/*add edge*/

		if (++m % 100000 == 0)
			DEBUG("Already %u lines read.\n", m);
    }
    stream.close();

	N = str2id.size();
	edges.resize(N);

	Nones = 0;
    for (int i=0; i<N; ++i)
		Nones += edges[i].size();

	long total_edges = (directed ? (long)N*(N-1) : (long)N*(N-1)/2);

	Nzeros = total_edges - Nones;

	ratio0 = (double)Nzeros / total_edges;
	ratio1 = (double)Nones / total_edges;

	id2str.resize(N);
	for (auto& m: str2id)
		id2str[m.second] = m.first;

	network.set_size(N, N);

	for (int i=0, n=0; i<N; ++i)
		for (auto iter=edges[i].begin(); iter!=edges[i].end(); ++iter)
		{
			int j = *iter;
			network(i, j) = 1;
		}

    INFO("Reading graph finished! [nodes:%d edges:%u]\n", N, Nones);

	return 0;
}

int Graph::check_id(const string& str)
{ 
    if (str2id.find(str) == str2id.end())
        str2id[str] = str2id.size();
    return str2id[str];
}

Graph::Heldout* Graph::create_heldouts(int* sizes[2], int num)
{
	INFO("Start assigning %d heldout sets!\n", num);
	long total_links = 0, total_nlinks = 0;

	for (int i = 0; i < num; i++)
	{
		total_nlinks += sizes[0][i];
		total_links += sizes[1][i];
	}

	// check
	int ratio = 10;
	if (total_links >= Nones/ratio || total_nlinks >= Nzeros/ratio)
	{
		ERROR("Too many links or nonlinks!\n");
		return NULL;
	}

	Heldout *heldouts = new Heldout[num];
	
	Heldout::pair_set pools[2];

	/* assigning non links */
    srand(time(0));
    while (total_nlinks-- > 0)
	{
        int source, dest;
        do
		{
        	source = rand()%N;
        	dest = rand()%N;
			if (!directed && source > dest)
				swap(source, dest);
        } while (source == dest
				|| pools[0].count(make_pair(source,dest)) >= 1 
				|| network(source, dest) != 0);

		pools[0].insert(make_pair(source,dest));
    }

	DEBUG("Finished assigning nonlinks.\n");

	set<unsigned> eindexes;
	while (eindexes.size() < total_links)
		eindexes.insert(rand() % Nones);

	auto nit = network.begin();
	unsigned pre = 0;
	for (set<unsigned>::iterator it=eindexes.begin(); it!=eindexes.end(); it++)
	{
		unsigned adv = *it-pre;
		pre = *it;
		while (adv--)
			++nit;
		int i = nit.row(), j = nit.col();
		pools[1].insert(make_pair(i, j));
	}
	for (const pair<int,int>& pr : pools[1])
	{
		network(pr.first, pr.second) = 0;
	}

	DEBUG("Finished assigning links.\n");

	for (int i=0; i<2; i++)
	{
		int bin = 0;
		for (const pair<int,int>& pr : pools[i])
		{
			if (heldouts[bin].pairs[i].size() == sizes[i][bin])
				bin++;
			heldouts[bin].pairs[i].insert(pr);
		}
	}

	for (int bin = 0; bin < num; bin++)
	{
		heldouts[bin].ratio1 = ratio1;
		heldouts[bin].ratio0 = ratio0;
	}

	INFO("Finished assigning %d heldout sets!\n", num);

    return heldouts;
}
