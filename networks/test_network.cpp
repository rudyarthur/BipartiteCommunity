#include <iostream>
#include "network_utils.h"
#include "network.h"

using namespace std; 

int main(int argc, char **argv){

	
	cout << ">> Empty network" << endl;
	Network G;
	G.print_basic();
	
	cout << ">> read network from file" << endl;
	print_header("../test_data/graph3.txt");
	G.read("../test_data/graph3.txt");
	G.print_basic();
	cout << ">> Adjacency Matrix" << endl;
	G.print_adjacency();
	
	cout << ">> Network as adjacency map" << endl;
	print_header("../test_data/graph3.txt");
	adjacency_map tmp = G.to_map();
	for(auto &i : tmp){
		cout << i.first << " :: ";
		for(auto &j: i.second){
			cout << "(node:" << j.first << ", weight:" << j.second << ") ";
		} cout << endl;
	}
	
	cout << ">> Network selfloops" << endl;
	for(auto &n : G.get_nbrs){
		cout << "G.selfloop(" << n.first << ") = " << G.selfloop(n.first) << endl;
	}

	cout << ">> Network join" << endl;
	print_header("../test_data/graph3.txt");
	cout << "Join 0 -> 2" << endl;
	G.join(2,0);
	cout << "Neighbours of 2 :: ";
	for(auto &j: G.get_nbrs[2]){
		cout << "(node:" << j.first << ", weight:" << j.second << ") ";
	} cout << endl;
	cout << "Join 4 -> 5" << endl;
	G.join(5,4);
	cout << "Neighbours of 5 :: ";
	for(auto &j: G.get_nbrs[5]){
		cout << "(node:" << j.first << ", weight:" << j.second << ") ";
	} cout << endl;
	G.print_connections();

	cout << ">> Network reindex" << endl;
	G.reindex();
	G.print_indexmap();
	cout << ">> Print reindexed network with original indices" << endl;
	G.print_connections(true);
	

	
	return 0;
}
