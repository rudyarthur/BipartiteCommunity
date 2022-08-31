#include <iostream>
#include "network_utils.h"
#include "weighted_network.h"

using namespace std; 

int main(int argc, char **argv){

	
	cout << ">> Empty network" << endl;
	WeightedNetwork G;
	G.print_basic();
	
	cout << ">> read network from file" << endl;
	print_header("../test_data/graph3.txt");
	G.read("../test_data/graph3.txt");
	G.print_basic();

	cout << ">> construct an induced network" << endl;
	print_header("../test_data/graph3.txt");
	unordered_map<int,int> labels{ {0,1},{1,1},{2,1},{3,2},{4,2},{5,2} };
	for(auto &i : labels){ cout << "node " << i.first << " has label " << i.second << endl; } 
	WeightedNetwork GI = G.induced_graph(labels);
	GI.print_connections();
	
	cout << ">> reindex the induced network" << endl;
	GI.reindex();
	GI.print_indexmap();
	GI.print_connections();
	
	return 0;
}
