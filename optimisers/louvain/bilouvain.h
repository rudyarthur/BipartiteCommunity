#ifndef BILOUVAIN_H
#define BILOUVAIN_H


#include "../optimiser.h"
#include "aggregate.h"

using namespace std; 

class BiLouvain : public Optimiser<BarberModularity> {
	public:
		bool max_agg;

		BiLouvain(BarberModularity& Q);		
		double optimise(); 
};


BiLouvain::BiLouvain(BarberModularity& _Q) : Optimiser<BarberModularity>(_Q) {
		optimiser_type = "bilouvain";
		max_agg = true;		
}


double BiLouvain::optimise() {
	
	double dq = 0;
	int level = 0;
	vector< unordered_map< int, int > > community_maps;
	int side = (2*Q.G.B > Q.G.num_nodes) ? 0 : 1; //do the bigger side first
	BipartiteNetwork Binit( Q.G );

	do{

		//find new community labels	
		dq = greedy_update(); 

		//save community assignment
		community_maps.push_back( Q.node_to_comm );
		
		//relabel to avoid annoying traps, what a nightmare!
		unordered_map<int, int> lab_map;
		unordered_map<int, int> labs;
		for(auto &i : Q.comm_to_node){ 
			for(auto &j : i.second){
				if( Q.G.colour[j] == side){lab_map[i.first] = j;}
			}
		}
		for(auto &i  : Q.node_to_comm){ 
			Q.node_to_comm[i.first] = lab_map[i.second]; 
			if( Q.G.colour[i.first] == side){ labs[i.first] = Q.node_to_comm[i.first];}
		}
		fixupQ();
		
		BipartiteNetwork BI = Q.G.induced_graph( labs , side ); 
		
		//reset quality function with GI
		Q.setup_network(BI);
		setup_sizes();
		side = 1 - side;
		
		++level;
		
		
	} while(dq > eps);

	assign_labels(community_maps);	
	
	{
		Q.G = Binit;
		for(auto &i: labels){ Q.node_to_comm[i.first] = labels[ i.first ]; }
		fixupQ();
		
		Aggregate< BarberModularity > A( Q, max_agg );
		A.optimise();

		labels.clear();
		for(auto &i: A.labels){ labels[ i.first ] = i.second; }
	}
	
	
	return 0;
}

#endif
