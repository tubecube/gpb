#include "graph.hpp"

int Graph::read_from_file(const string& filename, bool directed)
{
    ifstream stream(filename);
    if (stream.is_open() == false)
	{
       ERROR("Graph: fail opening file %s\n", filename.c_str());
	   return 2;
	}

    this->directed = directed;

    if (directed)
       INFO("Graph: reading from directed graph %s\n", filename.c_str());
    else
       INFO("Graph: reading from undirected graph %s\n", filename.c_str());

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

		check_edge(sid, did);
		edges[sid].insert(did);
		/*add edge*/

		if (++m % 100000 == 0)
			DEBUG("Graph: already %u lines read\n", m);
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

    INFO("Graph: reading graph finished nodes:%d edges:%u\n", N, Nones);

	return 0;
}

int Graph::check_id(const string& str)
{ 
    if (str2id.find(str) == str2id.end())
        str2id[str] = str2id.size();
    return str2id[str];
}

int Graph::pop_heldout()
{
	if (!heldouts.empty())
	{
		heldouts.pop_back();
		return 0;
	} return -1;
}

const Graph::Heldout* Graph::push_heldout(int N0, int N1)
{
	DEBUG("Graph: creating heldout with %d links and %d nonlinks\n", N1, N0);
	/*check*/
	unsigned total_nonlinks = N0;
	unsigned total_links = N1;
	for (const Heldout& heldout : heldouts)
	{
		total_nonlinks += heldout.pairs[0].size();
		total_links += heldout.pairs[1].size();
	}

	if (total_nonlinks >= Nzeros/2 || total_links >= Nones/2)
	{
		ERROR("Graph: links or nonlinks exceed half, not allowed!\n");
		return NULL;
	}

	heldouts.push_back(Heldout());
	Heldout& current = heldouts.back(); 

    srand(time(0));

	/* assigning nonlinks */
    while (N0-- > 0)
	{
        int source, dest;
        do
		{
        	source = rand()%N;
        	dest = rand()%N;
			check_edge(source, dest);
        } while (source == dest || network(source, dest) != 0 || check_in_heldouts(source, dest, false));
		current.pairs[0].insert(make_pair(source, dest));
    }

	/* assigning links */
	int left = N1 - current.pairs[1].size();
	while (left)
	{
		set<unsigned> eindexes;
		while (eindexes.size() < left)
			eindexes.insert(rand() % Nones);

		auto nit = network.begin();
		unsigned pre = 0;
		for (set<unsigned>::iterator it=eindexes.begin(); it!=eindexes.end(); it++)
		{
			unsigned adv = *it-pre;
			pre = *it;
			while (adv--)
				++nit;
			int source = nit.row();
			int dest = nit.col();
			if (!check_in_heldouts(source, dest, true))
				current.pairs[1].insert(make_pair(source, dest));
		}
		left = N1 - current.pairs[1].size();
	}

	return &current;
}

bool Graph::check_in_heldouts(int source, int dest, bool link) const
{
	int idx = link ? 1 : 0;
	check_edge(source, dest);
	pair<int,int> pr = make_pair(source, dest);
	for (const Heldout& heldout : heldouts)
		if (heldout.pairs[idx].count(pr) >= 1)
			return true;
	return false;
}
