/*
 * Harvester.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: Sergey Ivanov
 */
// Find seeds using Harvester
// Input: <filename of the weighted graph> <number of vertices> <folder with PWs> <number of seeds> <filename of the seeds>
// Example: datasets/hep15233.txt 15233 PW/MC/ 20 seeds/seeds.txt
// currently uses a series of folders PW/Flickr%s.txt and write to corresponding seeds files.

#include <iostream>
#include <fstream>					// for reading files
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <sys/time.h>
#include "dirent.h" 				// for reading files from directory
#include <queue> 					// std::priority_queue
#include <functional>				// std::greater
#include <sys/stat.h>				// to check for files
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
typedef typename boost::graph_traits<DeterministicGraph>::vertex_iterator vertex_iterator;
typedef typename boost::graph_traits<DeterministicGraph>::vertex_descriptor Vertex;
typedef typename graph_traits<WeightedGraph>::edge_descriptor Edge;
typedef typename property_map < WeightedGraph, edge_weight_t >::type myWeight;

struct CC {
	int idx;
	std::vector<int> nodes;
};

//			comparator
class CC_greater {
public:
	bool operator()(CC& cc1, CC& cc2)
	{
		return (cc1.nodes.size() > cc2.nodes.size());
	}
};
typedef std::priority_queue <CC, vector<CC>, CC_greater> CC_pq;

std::vector<Edge> edgeList;
double P=0;
int E=0;

// print execution time
double print_time(struct timeval &start, struct timeval &end){
	long double usec;

	usec = (end.tv_sec*1000000 + end.tv_usec) - (start.tv_sec*1000000 + start.tv_usec);
	return usec/1000000.0;
}

// read Graph with probabilities from file
void buildList (WeightedGraph& G, std::string input){

	typedef typename WeightedGraph::edge_property_type Weight;
	typename graph_traits<WeightedGraph>::edge_iterator ei, eiend;

	std::map<int, int> m;

	std::ifstream infile(input.c_str());
	if (infile==NULL){
		std::cout << "Unable to open the input file\n";
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
			std::cout << "Unable to insert edge\n";
		}
		edgeList.push_back(e);
		P+=p;
		edge_count++;
	}
	std::cout <<"E = "<<edge_count<< " P = " << P << " V = "<<mapped<<std::endl;
	E=edge_count;
}

// read Possible World from file
bool buildPW(DeterministicGraph& G, int V, string filename){
	typedef typename graph_traits<DeterministicGraph>::edge_descriptor Edge;
	typename graph_traits<DeterministicGraph>::edge_iterator ei, eiend;


	std::ifstream infile(filename.c_str());
	if (infile==NULL){
		std::cout << "Unable to open the input file "<<filename<<"\n";
	}
	Edge e;
	bool inserted;
	std::string line;
	int i =0;
	int count=0;
	int edge_count = 0;

	/* While there is still a line. */
	while(getline(infile, line)) {
		tokenizer<> tok(line);
		int neighbors = 0;
		for(tokenizer<>::iterator beg=tok.begin(); beg!=tok.end();++beg){
			neighbors++;
			int v=lexical_cast<int>(*beg);
			if (i<v){
				edge_count++;
				tie(e,inserted)=add_edge(i, v, G);
				if (!inserted) {
					std::cout << "Unable to insert edge\n";
				}
				else count++;
			}
		}
		// add isolated vertex to the graph
		// quick and dirty solution
		// for more advanced look at: http://stackoverflow.com/a/27648786/2069858
//		if (neighbors == 0) {
//			tie(e,inserted)=add_edge(i, i, G);
//			if (!inserted) {
//				std::cout << "Unable to insert edge\n";
//			}
//			remove_edge(i, i, G);
//		}
		i++;
	}
//	cout<< "Edges loaded= "<<count<<"\n";
	infile.close();
	cout << "n: " << i << " m: " << edge_count << endl;
	if (i == V)
		return true;
	else
		return false;
}

// read Possible World from file as adjacency list
bool readPW(vector<int> pw[], int V, string filename){
	ifstream infile(filename.c_str());
	if (infile==NULL){
		cout << "Unable to open the input file "<<filename<<"\n";
	}
	string line;
	int i = 0;
	int edge_count = 0;

	/* While there is still a line. */
	while(getline(infile, line)) {
		tokenizer<> tok(line);
		for(tokenizer<>::iterator beg=tok.begin(); beg!=tok.end();++beg){
			int v = lexical_cast<int>(*beg);
			if (i < v) {
				edge_count++;
				pw[i].push_back(v);
				pw[v].push_back(i);
			}
		}
		i++;
	}
	infile.close();
//	cout << "n: " << i << " m: " << edge_count << endl;
	if (i == V)
		return true;
	else
		return false;
}

// read Possible World from file as edge list
bool readPW2(vector<vector<int> >& pw, int V, string filename){
	ifstream infile(filename.c_str());
	if (infile==NULL){
		cout << "Unable to open the input file "<<filename<<"\n";
	}

	string line;
	int edge_count = 0;
	int u, v;
	double p;
	bool added = false;

	// makes an assumption that edge list includes an edge only once (eg., (1 3) or (3 1))
	while (infile >> u >> v >> p) {
//		cout << u << " " << v << endl;
		edge_count++;
		pw[u].push_back(v);
		pw[v].push_back(u);
		added = true;
	}
//	int count = 0;
//	for (int i=0; i < V; ++i) {
//		if (pw[i].size() > 250) {
//			count++;
//			cout << i << ": ";
//			for (int j = 0; j < pw[i].size(); ++j)
//				cout << pw[i][j] << " ";
//			cout << endl;
//		}
//	}
//	cout << count << endl;

	infile.close();
//	cout << "n: " << i << " m: " << edge_count << endl;
	if (added)
		return true;
	else
		return false;
}

class compare {
public:
	bool operator()(vector<int> lhs, vector<int> rhs) {
		return lhs.size() > rhs.size();
	}
};

typedef priority_queue<vector<int>, vector<vector<int> >, compare> PQ;

// find top-k CC using boost graph representation
PQ find_top_CC(DeterministicGraph& PW, int k) {
	typename graph_traits <DeterministicGraph>::out_edge_iterator out1, out2;

	PQ pq; // priority queue for top-k CC
	int V = num_vertices(PW);
	bool *explored = new bool[V+1];
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
				Vertex u = vertex(CC[ix], PW);
				for (tie(out1, out2) = out_edges(u, PW); out1 != out2; ++out1) {
					Vertex v = target(*out1, PW);
					if (!explored[v]) {
						CC.push_back(v);
						explored[v] = true;
					}
				}
				ix++;
			}
//			if (CC.size() > 5)
//				cout << cc_number << ": " << CC.size() << endl;
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
//	cout << "Total CCs: " << cc_number << endl;
//	cout << "Unit CCs: " << unit_sized << endl;
	return pq;
}

// find top-k CCs using vectors graph representation
PQ find_top_CC2(vector<vector<int> >& PW, int V, int k) {

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
//			if (CC.size() > 5)
//				cout << cc_number << ": " << CC.size() << endl;
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
//	cout << "Total CCs: " << cc_number << endl;
//	cout << "Unit CCs: " << unit_sized << endl;
	return pq;
}

int main(int argc, char* argv[]) {
	struct timeval ex_start, ex_finish;
	gettimeofday(&ex_start, NULL);

	// read parameters from command-line
	const std::string dataset = argv[1]; // filename of the dataset
	const long int V = atoi(argv[2]); // number of nodes
	const std::string pw_dir = argv[3]; // folders containing PWs
	const int k = atoi(argv[4]);         // number of seeds
	const std::string outfilename = argv[5]; // filename of the seeds

	cout << "Started processing..." <<endl;

	cout << "Graph: " << dataset << " k: " << k << endl;

	// read graph from the file
	WeightedGraph G(V);

	struct timeval t_start, t_finish;
	gettimeofday(&t_start, NULL);
	buildList(G, dataset);
	gettimeofday(&t_finish, NULL);
	cout << "* Read graph: " << print_time(t_start, t_finish) << "sec." << endl;
//	return 1;

	// update scores
	struct timeval update_start, update_finish;
	gettimeofday(&update_start, NULL);

//	for sparsification graphs
//	std::ostringstream s;
//	s << "PW/G" << K << "/";
//	std::string pw_folder(s.str());
//	cout << pw_folder << endl;

	// For every filename in a PW directory create graph PW and update scores accordingly
	// solution found here: http://www.dreamincode.net/forums/topic/59943-accessing-directories-in-cc-part-i/
	DIR *pdir = NULL;
//	pdir = opendir(pw_folder.c_str());
    pdir = opendir(pw_dir.c_str());
	struct dirent *pent = NULL;
	struct stat filestat;

	if (pdir == NULL) // if pdir wasn't initialised correctly
	{ // print an error message and exit the program
		cout << "\nERROR! pdir could not be initialised correctly";
		exit (3);
	} // end if

	std::map <int, double> scores;
	for (int i = 0; i < V; ++i) {
		scores[i] = 0;
	}
	int PWs = 0;
	double iter_time = 0;
	double pw_time = 0;
	double update_time = 0;
	double topcc_time = 0;

	// while there is still something in the directory to list
	while ( (pent = readdir(pdir)) != NULL ) {

		// if pent has not been initialised correctly
		if (pent == NULL) { // print an error message, and exit the program
			cout << "\nERROR! pent could not be initialised correctly";
			exit(3);
		}

		  string filepath = pw_dir + "/" + pent->d_name;
//		 string filepath = pw_folder + "/" + pent->d_name;

		 // check that files are text ones
		 // http://www.cplusplus.com/forum/beginner/10292/
		 if (stat(filepath.c_str(), &filestat)) continue;
		 if (S_ISDIR( filestat.st_mode ))         continue;


		struct timeval t_start, t_finish;
		gettimeofday(&t_start, NULL);
//        DeterministicGraph PW[V];
//        cout << "\n" << "New PW: ";
//        cout << pw_dir + pent->d_name << endl;
//        buildPW(*PW, pw_dir + pent->d_name);
//        vector<int>* PW = new vector<int>[V];
		vector<vector<int> > PW(V);
		 bool built = readPW2(PW, V, pw_dir + pent->d_name);
//		bool built = readPW2(PW, V, pw_folder + pent->d_name);

//        DeterministicGraph PW[V];
//        bool built = buildPW(*PW, V, pw_dir + pent->d_name);

//		vector<int> PW[V];
//		bool built = readPW(PW, V, "PW/MC/MC1_Sergei.txt");
//		DeterministicGraph PW[V];
//		bool built = buildPW(*PW, V, "PW/MC/MC1_Sergei.txt");


//        bool built = readPW(PW, V, "PW/MC/MC1_Sergei.txt");
//        cout << "nodes: " << num_vertices(*PW) << " edges: " << num_edges(*PW) << endl;
		gettimeofday(&t_finish, NULL);
		pw_time += print_time(t_start, t_finish);
//        cout << "* Read PW: " << print_time(t_start, t_finish) << "sec." << endl;

		if (built) {
			PWs++;
//        	cout << PWs << ": " << pent->d_name << " ";
			struct timeval iter_start, iter_finish;
			gettimeofday(&iter_start, NULL);

			// find scores
			struct timeval topcc_start, topcc_finish;
			gettimeofday(&topcc_start, NULL);
			PQ pq;
//			pq = find_top_CC(*PW, k);
			pq = find_top_CC2(PW, V, k);
			gettimeofday(&topcc_finish, NULL);
			topcc_time += print_time(topcc_start, topcc_finish);

			// update scores
			struct timeval assign_start, assign_finish;
			gettimeofday(&assign_start, NULL);

			for (int j=0; j<k; ++j) {
				vector<int> top = pq.top();
				pq.pop();
				int top_size = top.size();
				for (int ix = 0; ix < top_size; ++ix) {
					scores[top[ix]] += 1./top_size;
//					cout << "node " << top[ix] << " score " << scores[top[ix]] << endl;
				}
//				cout << " size " << top.size() << endl;
			}

//			for (int j=0; j<k; ++j) {
//				CC top = pq.top();
//				pq.pop();
//				int top_size = top.nodes.size();
//				for (int ix = 0; ix < top_size; ++ix) {
//					scores[top.nodes[ix]] += 1./top_size;
////					cout << "node " << top.nodes[ix] << " score " << scores[top.nodes[ix]] << endl;
//				}
////				cout << "cc " << top.idx << " size " << top.nodes.size() << endl;
//			}
			gettimeofday(&assign_finish, NULL);
			update_time += print_time(assign_start, assign_finish);

//			break;
			gettimeofday(&iter_finish, NULL);
			iter_time += print_time(iter_start, iter_finish);
//			cout << print_time(iter_start, iter_finish) << endl;
//			if (PWs == 3)
//				break;
		}
	}
	gettimeofday(&update_finish, NULL);
	cout << "* Accumulation phase: " << print_time(update_start, update_finish) << "sec." << endl;
//    cout << "  Average time per PW: " << iter_time/PWs << "sec." << endl;
	cout << "  ** Read PW: " << pw_time << "sec." << endl;
//    cout << "  ** Discover CC per PW: " << cc_time << "sec." << endl;
//    cout << "  ** Maintain PQ per PW: " << pq_time << "sec." << endl;
	cout << "  ** Find top-k CCs: " << topcc_time << "sec." << endl;
	cout << "  ** Assign scores: " << update_time << "sec." << endl;

	// select nodes
	struct timeval select_start, select_finish;
	gettimeofday(&select_start, NULL);

	std::set<int> S;
	std::map <int, bool> selected;
	for (int i = 0; i < V; ++i) {
		selected[i] = false;
	}

	typename graph_traits <WeightedGraph>::out_edge_iterator out1, out2;
	myWeight weight;

	for (int ix = 0; ix < k; ++ix) {
		// find max element
		int max_node;
		double max_score = 0;
		for (std::map<int, double>::iterator it = scores.begin(); it != scores.end(); ++it) {
//			cout << "inside for loop" << endl;
			if (it->second > max_score) {
				max_score = it->second;
				max_node = it->first;
			}
		}
//		cout << endl;
//		cout << "max_node: " << max_node << " max_score: " << max_score << endl;
		S.insert(max_node);
		selected[max_node] = true;
		scores.erase(max_node);

		// penalize scores of neighbors
		Vertex u = vertex(max_node, G);
		for (tie(out1, out2) = out_edges(u, G); out1 != out2; ++out1) {
			Vertex v = target(*out1, G);
			if (!selected[v]) {
				double p = get(weight, *out1);
				scores[v] *= (1 - p);
			}
//			cout << u << "-->" << v << endl;
		}
	}
	gettimeofday(&select_finish, NULL);
	cout << "* Penalization phase: " << print_time(select_start, select_finish) << "sec." << endl;

	for (std::set<int>::iterator it = S.begin(); it != S.end(); ++it) {
		std::cout << *it << ' ';
	}
	cout << endl;

	// write seeds to file
	struct timeval write_start, write_finish;
	gettimeofday(&write_start, NULL);

//	for sparsified graphs
//	std::ostringstream seeds_name;
//	seeds_name << "seeds/seeds" << K << ".txt";
//	std::string seeds_filename(seeds_name.str());

	std::ofstream myfile;
	myfile.open(outfilename.c_str());
//	myfile.open(seeds_filename.c_str());
//	cout << seeds_filename;
//	myfile << k << " ";
	for (std::set<int>::iterator it = S.begin(); it != S.end(); ++it) {
		myfile << *it << ' ';
	}
	myfile.close();

	gettimeofday(&write_finish, NULL);
	cout << "* Wrote seeds to file: " << print_time(write_start, write_finish) << "sec." << endl;

	closedir(pdir);
	gettimeofday(&ex_finish, NULL);

	cout << endl;
	cout << "Total number of PWs: " << PWs << endl;

	cout << "* Execution time of Harvester: " << print_time(ex_start, ex_finish) << " sec." << endl;

	cout << "\nEnd of Harvester" << endl;
}
