#include "graph.hpp"

void
Graph::read_from_file(const string& filename, const string& ground, bool directed, bool heldout)
{
    this->directed = directed;
    if (directed)
        cout << "reading directed network\n";
    else
        cout << "reading undirected network\n";

    ifstream stream(filename);
    if (stream.is_open() == false)
    {
        cerr << " fail opening file !\n";
        exit(-1);
    }
    for (string line; getline(stream, line); )
    {
        istringstream fields(line);
        string first, second;
        fields >> first >> second;
        add_edge(first, second);
    }
    stream.close();

    // 有可能有重复边 
    M = 0;
    for (int i=0; i<N; ++i)
        M += network[i].size();

    cout << M << " edges" << endl;
    cout << N << " nodes" << endl;

    // create ground truth communities
    create_ground_truth(ground);

    if (heldout) create_heldout(M*0.05, M*0.05);

    if (directed) network2.resize(N);
    for (int i=0; i<N; ++i)
        for (auto it=network[i].begin(); it!=network[i].end(); ++it)
            network2[*it].insert(i);

    // 清空str2id, 创建id2str
    id2str.resize(N);
    for (auto it=str2id.begin(); it!=str2id.end(); ++it)
        id2str[it->second] = it->first;
    str2id.clear();

}

void Graph::create_ground_truth(const string& ground)
{
    ifstream stream(ground);
    if (stream.is_open())
    {
        for (string line; getline(stream, line); )
        {
            set<int> s;
            istringstream fields(line);
            for (string field; fields >> field; )
            {
                if (field.length() == 0)
                    cerr << "Warning: two consecutive tabs, or tab at the start of a line. Ignoring empty fields like this" << endl;
                else
                {
                    if (str2id.find(field) == str2id.end())
                        cerr << "Warning: " << field << " not exist in input graph" << endl;
                        // continue;
                    else {
                        s.insert(str2id[field]);
                        ground_set.insert(str2id[field]);
                    }
                }
            }
            if (s.size() == 0)
                cerr << "Warning: ignoring empty sets in file: " << ground << endl;
            else
                ground_truth.push_back(s);
        }
    }
}

void
Graph::add_edge(string start, string end)
{
    if (start.size() == 0 || end.size() == 0)
        return;
    if (start[0] == '%') {
        cout << "comment";
        return;
    }
    int startid = get_id(start);
    int endid = get_id(end);
    if (network.size() < N)
        network.resize(N);

    if (directed) {
        network[startid].insert(endid);
    } else {
        if (startid < endid)
            network[startid].insert(endid);
        else
            network[endid].insert(startid);
    }
    M++;
    if (M % 10000 == 0)
        cout << M << " edges" << flush << '\r';
}

int
Graph::get_id(const string& str)
{ 
    int id;
    auto got = str2id.find(str);
    if (got == str2id.end())
    {
        id = str2id.size();
        str2id[str] = id;
        N++;
    }
    else
        id = got->second;
    return id;
}

map<pair<int,int>,bool> &
Graph::create_heldout(size_t size1, size_t size2)
{
    cout << "creating heldout set <" << size1 << "> <" << size2 << ">" <<endl;

    this->heldout.clear();

    srand(time(0));
    while(size2 > 0) {
        int i,j;
        do {
        i = rand()%N;
        j = rand()%N;
        } while (i==j || network[i].find(j) != network[i].end());
        auto ret = heldout.insert(pair<pair<int,int>,bool>(pair<int,int>(i,j), false));
        if (ret.second==true) size2--;
    }
    set<int> ran;
    while (size1 > 0) {
        int i = rand()%M;
        auto ret = ran.insert(i);
        if (ret.second==true) size1--;
    }

    auto it = ran.begin();
    auto it_end = ran.end();

    int ni = 0;
    for (int i=0; i<N; ++i) {
        for (auto itr=network[i].begin(); itr!=network[i].end(); ++itr) {
            if (ni == *it) {
                ++it;
                heldout.insert(pair<pair<int,int>,bool>(pair<int,int>(i,*itr), true));
            }
            if (it==it_end) break;
            ni++;
        }
        if (it==it_end) break;
    }
    // erase node pair from network which in heldout
    for (auto itr=heldout.begin(); itr!=heldout.end(); ++itr)
        if (itr->second == true) {
            int i = itr->first.first;
            int j = itr->first.second;
            network[i].erase(j);
            M--;
        }
    return heldout;
}

void Graph::save_community(vector<set<int>>& comm) {
    ofstream ofs("communities.txt");
    for (int k=0; k<comm.size(); ++k) {
        for (auto itr = comm[k].begin(); itr != comm[k].end(); ++itr)
            ofs << id2str[*itr] << " ";
        ofs << endl; 
    }
    ofs.close();
}

void Graph::cut(vector<set<int>>& group) {
    if (ground_set.size() == 0)
        return;
    for (int i=0; i<group.size(); ++i) {
        auto end = group[i].end();
        auto it = group[i].begin();
        while (it != end) {
            if (ground_set.find(*it) == ground_set.end())
                it = group[i].erase(it);
            else
                it++;
        }
    }
}
