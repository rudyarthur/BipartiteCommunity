#ifndef WEIGHTEDNETWORK_H
#define WEIGHTEDNETWORK_H

#include "network.h"
#include <unordered_map>
#include <unordered_set>

using namespace std; 

///An undirected weighted network
class WeightedNetwork : public Network {
	public:
				
		WeightedNetwork():Network(){
			network_type="weighted_network";
		}		
		
		WeightedNetwork(const char* filename);
		WeightedNetwork(adjacency_map &lk);
		WeightedNetwork(const Network &G );
	
		WeightedNetwork induced_graph(unordered_map<int, int> &label, bool unweighted);
};



///Use provided node labels to compute induced graph
/**
 * labels are the node labels!
 * unweighted = True computes unweighted adjacency, all edges have weight 1
 * I_ab = sum_{i, l(i) = a} sum_{j, l(j) = b} A_{ij} = sum_{ij} A_{ij} delta( l(i), a ) delta( l(j), b )  
**/
WeightedNetwork WeightedNetwork::induced_graph(unordered_map<int, int> &label, bool unweighted=false){ 
	if( (unsigned)num_nodes != label.size()){
		cerr << "Trying to induce a " << num_nodes << " network with " << label.size() << " labels" << endl;
		exit(1);
	}
	adjacency_map mp;
	for(auto &i : label){
		int a = i.second; //a is the label of i
		for(auto &j : get_nbrs[i.first]){
			int b = label[ j.first ];
			if(mp.find(a) == mp.end()){ mp[a][b] = (unweighted) ? 1 : j.second; }
			else{
				if(mp[a].find(b) == mp[a].end()){ mp[a][b] = (unweighted) ? 1 : j.second;  }
				else{ mp[a][b] += (unweighted) ? 0 : j.second; }
			}
		}
	}
	return WeightedNetwork(mp);
	
}


WeightedNetwork::WeightedNetwork(const char* filename){
	network_type="weighted_network";	
	read(filename);	
}

WeightedNetwork::WeightedNetwork(adjacency_map &mp) : Network(mp){
	network_type="weighted_network";	
}


WeightedNetwork::WeightedNetwork(const Network &G) : Network(G){
	network_type="weighted_network";	
}




#endif
