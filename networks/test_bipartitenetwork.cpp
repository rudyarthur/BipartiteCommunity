#include <iostream>
#include "network_utils.h"
#include "weighted_network.h"
#include "bipartite_network.h"

using namespace std; 

int main(int argc, char **argv){

	
	cout << ">> Empty network" << endl;
	BipartiteNetwork G;
	G.print_basic();
	
	cout << ">> read network from file" << endl;
	print_header("../test_data/bipartite1.txt");
	G.read("../test_data/bipartite1.txt");
	G.print_basic();
	

	
	cout << ">> construct an induced network" << endl;
	print_header("../test_data/bipartite1.txt");
	unordered_map<int,int> labels{ {0,0},{1,0},{2,0},{3,1},{4,2},{5,2},{6,0},{7,0},{8,1},{9,2} };
	for(auto &i : labels){ cout << "node " << i.first << " has label " << i.second << endl; } 
	BipartiteNetwork GI = G.induced_graph(labels, false, 0);
	GI.print_basic();
	
	cout << ">> reindex the induced network" << endl;
	GI.reindex();
	GI.print_indexmap();
	GI.print_connections();
	
	return 0;
}
