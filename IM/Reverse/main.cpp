#include <iostream>
#include <fstream>                  // for reading files
#include <algorithm>                 // for std::for_each
#include "dirent.h"                 // for reading files from directory
#include <queue>                    // std::priority_queue
#include <map>
#include <iomanip>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

//using namespace boost;
using namespace std;

struct CC {
    int idx;
    vector<int> nodes;
};

//          comparator
class CC_greater {
public:
    bool operator()(CC& cc1, CC& cc2)
    {
        return (cc1.nodes.size() > cc2.nodes.size());
    }
};

class compare {
public:
    bool operator()(vector<int> lhs, vector<int> rhs) {
        return lhs.size() < rhs.size();
    }
};

typedef priority_queue<vector<int>, vector<vector<int> >, compare> PQ;

// find top CCs that reach T nodes
vector<vector<int> > find_top_CC(vector<vector<int> >& PW, int V, int T) {

    PQ pq; // priority queue for top-k CC
    map<int, bool> explored;
    for (int i = 0; i < V; ++i)
        explored[i] = false;
    int cc_number = 0;

    for (int i = 0; i < V; ++i) {
        if (!explored[i]) {
            // perform BFS for vertex i
            vector<int> CC;
            CC.push_back(i);
            explored[i] = true;
            vector<int>::size_type ix = 0;
            while (ix < CC.size()) {
                int pw_size = PW[CC[ix]].size();
                for (int j = 0; j < pw_size; ++j) {
                    int v = PW[CC[ix]][j];
                    if (!explored[v]) {
                        CC.push_back(v);
                        explored[v] = true;
                    }
                }
                ix++;
            }
            cc_number++;

            pq.push(CC);
        }
    }
    int coverage = 0;
    vector<int> CC;
    vector<vector<int> > topT_CC;
    while (coverage < T) {
        CC = pq.top();
        coverage += CC.size();
        pq.pop();
        topT_CC.push_back(CC);
    }
    return topT_CC;
}

bool buildPW(vector<vector<int> >& PW, string dataset) {
    ifstream infile(dataset.c_str());
    if (infile==NULL){
        cout << "Unable to open the input file "<<dataset<<"\n";
    }

    int u, v;
    double p, r;
    // makes an assumption that edge list includes an edge only once (eg., (1 3) or (3 1))
    while (infile >> u >> v >> p) {
        r = ((double) rand() / (RAND_MAX));
        if (r < p) {
            PW[u].push_back(v);
            PW[v].push_back(u);
        }
    }
    infile.close();
    return true;
}

void readGraph(vector<vector<pair<int, double> > >& G, string dataset_f) {
    ifstream infile(dataset_f.c_str());
    if (infile==NULL){
        cout << "Unable to open the input file\n";
    }
    int u, v;
    double p;
    pair<int, double> endpoint1, endpoint2;
    vector<int> a;
    while (infile >> u >> v >> p){
        endpoint1 = make_pair(u, p);
        endpoint2 = make_pair(v, p);
        G[v].push_back(endpoint1);
        G[u].push_back(endpoint2);
    }
}

int accumulate_scores(map<int, double>& scores, string dataset, int T, int V, int R) {

    float amount = 0;
    for (int i=0; i<R; i++) {
        vector<vector<int> > PW(V);
        vector<vector<int> > topT_CC;
        bool built = buildPW(PW, dataset);
        if (built) {
            topT_CC = find_top_CC(PW, V, T);
            for (int j=0; j<topT_CC.size(); j++) {
                for (int ix=0; ix<topT_CC[j].size(); ++ix) {
                    int node = topT_CC[j][ix];
                    scores[node] += 1.0/topT_CC[j].size();
                }
            }
        }
        amount += (float)topT_CC.size()/(float)R;
    }
    int L = ceil(amount);

    return L;
}

void select_nodes(vector<int>& S, map<int, double>& scores, map<int, bool>& selected, int L, vector<vector<pair<int, double> > > G) {

    int argmax_v;
    double max_score;
    for (int i=0; i<L; ++i) {
        argmax_v = -1;
        max_score = 0;
        for (int j=0; j<scores.size(); j++) {
            if (!selected[j]) {
                if (scores[j] > max_score) {
                    max_score = scores[j];
                    argmax_v = j;
                }
            }
        }
        if (argmax_v == -1) {
            cout << "Failed to find next node..." << endl;
            throw invalid_argument("Failed to find next node...");
        }
        else {
            S.push_back(argmax_v);
            selected[argmax_v] = true;
            for (vector<pair<int, double> >::iterator it=G[argmax_v].begin(); it!=G[argmax_v].end(); ++it) {
                int u = it->first;
                double p = it->second;
                if (!selected[u]) {
                    scores[u] *= (1-p);
                }
            }
        }
    }
}

double calculateSpread(vector<vector<pair<int, double> > > G, vector<int> S, int V, int I) {
    double cascade_size = 0;
    int iter = 0;
    map<int, bool> activated;
    vector<int> T;
    while (iter < I) {
        // activate seeds
        for (int i = 0; i < V; ++i) {
            activated[i] = false;
        }
        for (vector<int>::iterator it = S.begin(); it != S.end(); ++it) {
            activated[*it] = true;
            T.push_back(*it);
        }
        // activate new nodes
        vector<int>::size_type ix = 0; // activated node index
        int u;
        while (ix < T.size()) {
            u = T[ix];
            for (vector<pair<int, double> >::iterator it = G[u].begin(); it != G[u].end(); ++it) {
                int v = it->first;
                double p = it->second;
                if (!activated[v]) {
                    double random = (double) rand()/RAND_MAX;
                    if (random <= p) {
                        activated[v] = true;
                        T.push_back(v);
                    }
                }
            }
            ix++;
        }
        cascade_size += T.size();
        iter++;
        T.clear();
    }
    double spread = (double)cascade_size/(double)I;
    return spread;
}

void add_more_nodes(vector<int>& S, int T, vector<vector<pair<int, double> > > G, map<int, double> scores, int L, int V, int I, map<int, bool>& selected) {
    // calculate spread of selected initial nodes
    double sigma = calculateSpread(G, S, V, I);
    int High=L;
    while (sigma < T) {
        select_nodes(S, scores, selected, High, G);
        High *= 2;
        sigma = calculateSpread(G, S, V, I);
    }
    int Low = floor(High/2);
    int k;
    vector<int> resultS;
    while (Low+1 < High) {
        k = Low + floor((High - Low)/2);
        resultS = vector<int> (&S[0], &S[k]);
        sigma = calculateSpread(G, resultS, V, I);
        if (sigma >= T)
            High = k;
        else
            Low = k;
    }
    if (sigma >=T)
        S = vector<int> (&S[0], &S[Low]);
    else
        S = vector<int> (&S[0], &S[High]);
    cout << "Final nodes: ";
    for (int i=0; i<S.size(); i++) {
        cout << S[i] << " ";
    }
    sigma = calculateSpread(G, resultS, V, I);
    cout << "Spread: " << sigma << " Wanted spread:" << T << endl;
}

int main(int argc, char* argv[]) {
    srand(time(NULL));

    // read parameters from command-line
    const string dataset = argv[1]; // filename of the dataset
    const long int V = atoi(argv[2]); // number of nodes
    const int R = atoi(argv[3]); // number of PWs
    const int T = atoi(argv[4]);         // threshold of active nodes
    const int I = atoi(argv[5]);   // number of MC calculations
    const string outfilename = argv[6]; // filename of the seeds

    cout << "Started processing..." <<endl;

    cout << "Graph: " << dataset << " T: " << T << endl;

    // read graph from the file
    cout << "Reading graph form the file..." << endl;
    vector<vector<pair<int, double> > > G(V);
    readGraph(G, dataset);

    //    accumulate scores
    cout << "Accumulating scores..." << endl;
    map<int, double> scores;
    int L = accumulate_scores(scores, dataset, T, V, R);

    //    select initial nodes
    cout << "Selecting initial nodes..." << endl;
    vector<int> S;
    map<int, bool> selected;
    select_nodes(S, scores, selected, L, G);

    //    add more nodes
    cout << "Adding the rest nodes..." << endl;
    add_more_nodes(S, T, G, scores, L, V, I, selected);

    //    writing to file seed sizes
    ofstream result_file(outfilename, ios_base::app);
    result_file << T << " " << S.size() << endl;
    result_file.close();

    cout << "Finished code" << endl;
    return 0;
}