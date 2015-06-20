#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

void Categories(string dataset_f, long int N, string prob_dataset) {

    ifstream infile(dataset_f.c_str());
    if (infile==NULL){
        cout << "Unable to open the input file\n";
    }

    map<long int, long int> nodes;
    vector<vector<long int> > graph(N);

    long int  u, v;
    // long int edge_count=0;
    long int count = 0;


    cout << "Start reading the file..." << endl;
    while (infile >> u >> v){
        if (!nodes.count(u)) {
            nodes[u] = count;
            count++;
            if (count==1000000 or count==2000000 or count==N-1)
                cout << count << endl;
        }
        if (!nodes.count(v)) {
            nodes[v] = count;
            count++;
            if (count==1000000 or count==2000000 or count==N-1)
                cout << count << endl;
        }
        graph[nodes[u]].push_back(v);
        graph[nodes[v]].push_back(u);
        // edge_count++;
    }
	// cout <<"E = "<<edge_count<< " V = "<< count <<endl;
    infile.close();

    cout << "Preparing probabilities..." << endl;
    vector<long int>  sizes(N);
    for (long int  i=0; i<graph.size(); i++) {
        sizes[i] = graph[i].size();
    }
    sort(sizes.begin(), sizes.end());
    long int  idx1 = N/4, idx2 = N/2, idx3 = 3*N/4;

    cout << "Writing probabilities to file..." << endl;
    // write Categories model to file
    ofstream ofile(prob_dataset.c_str());
    if (ofile==NULL){
        cout << "Unable to open the output file\n";
    }
    infile.open(dataset_f.c_str());
    double p;
    while (infile >> u >> v){
        if (sizes[nodes[u]] > sizes[idx3] or sizes[nodes[v]] > sizes[idx3])
            p = 0.8;
        else if (sizes[nodes[u]] > sizes[idx2] or sizes[nodes[v]] > sizes[idx2])
            p = 0.4;
        else if (sizes[nodes[u]] > sizes[idx1] or sizes[nodes[v]] > sizes[idx1])
            p = 0.2;
        else
            p = 0.1;
        ofile << nodes[u] << " " << nodes[v] << " " << p << endl;
    }
    ofile.close();
    infile.close();

}

int main(int argc, char* argv[]) {

    const string dataset = argv[1]; // filename of the dataset
    const long int V = atoi(argv[2]); // number of nodes
    const string output = argv[3]; // filename of the dataset with probabilities

    Categories(dataset, V, output);
    return 0;
}