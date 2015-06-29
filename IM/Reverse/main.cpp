#include <iostream>
#include <fstream>                  // for reading files
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <sys/time.h>
#include "dirent.h"                 // for reading files from directory
#include <queue>                    // std::priority_queue
#include <functional>               // std::greater
#include <sys/stat.h>               // to check for files
#include <string>
#include <sstream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/functional/hash.hpp>
#include "boost/lexical_cast.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
//#include "boost/filesystem/operations.hpp"
//#include "boost/filesystem/path.hpp"
//#include "boost/progress.hpp"

using namespace boost;
using namespace std;
typedef property < edge_weight_t, double >Weight;
typedef adjacency_list < vecS, vecS, undirectedS, no_property, Weight > WeightedGraph;
typedef adjacency_list < vecS, vecS, undirectedS> DeterministicGraph;
typedef typename graph_traits<DeterministicGraph>::vertex_iterator vertex_iterator;
typedef typename graph_traits<DeterministicGraph>::vertex_descriptor Vertex;
typedef typename graph_traits<WeightedGraph>::edge_descriptor Edge;
typedef typename property_map < WeightedGraph, edge_weight_t >::type myWeight;

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
typedef priority_queue <CC, vector<CC>, CC_greater> CC_pq;

class compare {
public:
    bool operator()(vector<int> lhs, vector<int> rhs) {
        return lhs.size() > rhs.size();
    }
};

typedef priority_queue<vector<int>, vector<vector<int> >, compare> PQ;

// find top-k CCs using vectors graph representation
PQ find_top_CC(vector<vector<int> >& PW, int V, int k) {

    PQ pq; // priority queue for top-k CC
    map<int, bool> explored;
    for (int i = 0; i < V; ++i)
    explored[i] = false;
    int cc_number = 0;
    int unit_sized = 0;

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
            //          if (CC.size() > 5)
            //              cout << cc_number << ": " << CC.size() << endl;
            if (CC.size() == 1) {
                unit_sized++;
            }
            cc_number++;

            // maintain CC priority queue
            int pq_size = pq.size();
            if (pq_size == k) {
                vector<int> top = pq.top();
                if (top.size() < CC.size()) {
                    pq.pop();
                    pq.push(CC);
                }
            }
            else
                pq.push(CC);
        }
    }
        //  cout << "Total CCs: " << cc_number << endl;
        //  cout << "Unit CCs: " << unit_sized << endl;
    return pq;
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

// read Graph with probabilities from file
void buildGraph (WeightedGraph& G, string input){

    typedef typename WeightedGraph::edge_property_type Weight;
    typename graph_traits<WeightedGraph>::edge_iterator ei, eiend;

    map<int, int> m;

    ifstream infile(input.c_str());
    if (infile==NULL){
        cout << "Unable to open the input file\n";
    }
    int u, v, u_, v_;
    double p;
    Edge e;
    bool inserted;
    int mapped=0;
    int edge_count=0;

    while (infile >> u >> v >> p){
        if (m.find(u)== m.end()){
            m[u]=mapped;
            mapped++;
        }
        if (m.find(v)==m.end()){
            m[v]=mapped;
            mapped++;
        }
        u_=m[u]; v_=m[v];
        tie(e,inserted)=add_edge(u_, v_, Weight(p), G);
        if (!inserted) {
            cout << "Unable to insert edge\n";
        }
        edge_count++;
    }
    cout <<"E = "<< edge_count << " V = "<< mapped << endl;
}

void accumulate_scores(string dataset, map<int, double>& scores, int& L, int T, int V, int R) {
    cout << "In the accumulate function" << endl;
}

int main(int argc, char* argv[]) {

    // read parameters from command-line
    const string dataset = argv[1]; // filename of the dataset
    const long int V = atoi(argv[2]); // number of nodes
    const int R = atoi(argv[3]); // number of PWs
    const int T = atoi(argv[4]);         // threshold of active nodes
    const string outfilename = argv[5]; // filename of the seeds

    cout << "Started processing..." <<endl;

    cout << "Graph: " << dataset << " T: " << T << endl;

    // read graph from the file
    WeightedGraph G(V);
    buildGraph(G, dataset);

    //    accumulate scores
    map<int, double> scores;
    int L;
    accumulate_scores(dataset, scores, L, T, V, R);

    //    select initial nodes
    vector<int> S;
    select_initial_nodes(S, scores, L, G);

    //    add more nodes
    add_more_nodes(S, T, G, scores);

    cout << "Finished code" << endl;
    return 0;
}