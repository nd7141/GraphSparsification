/*
 * main.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: sivanov
 */
// usage: ./main dataset R PW/,
// dataset is the filename of dataset
// R is the number of PWs to generate
// PW/ is the folder to save them.

#include <iostream>
#include <fstream> // for reading files
#include <cstdlib> // for atoi, rand
#include <map>
#include <tuple> // for tie
#include <sys/time.h> // for gettimeofday
#include "dirent.h" // for reading files from directory
#include <stdio.h> // printf()
#include <stdlib.h> // exit()
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstring>

using namespace std;

int main(int argc, char* argv[]) {

    const std::string dataset = argv[1]; // filename of the dataset
	const long int R = atoi(argv[2]); // number of PWs to generate
	const std::string pw_dir = argv[3]; // folder to write PW

	ifstream infile(dataset.c_str());
	if (infile==NULL){
		cout << "Unable to open the input file "<<dataset<<"\n";
	}
	string line;
	int edge_count = 0;
	int u, v;
	double p, r;
	string PWs[R];

	// makes an assumption that edge list includes an edge only once (eg., (1 3) or (3 1))
	while (infile >> u >> v >> p) {
		edge_count++;
		for (int i=0; i<R; i++) {
			r = ((double) rand() / (RAND_MAX));
			if (r < p) {
				ostringstream prob;
				prob << p;
				PWs[i] += to_string(u) + " " + to_string(v) + " " + prob.str() + "\n";
			}
		}
	}
	infile.close();

//	write PWs to files
	ofstream myfile;
	string pw_filename;
	for (int i=0; i<R; i++) {
		pw_filename = pw_dir + to_string(i+1) + ".txt";
		myfile.open(pw_filename);
		myfile << PWs[i];
		myfile.close();
	}

	return 0;
}




