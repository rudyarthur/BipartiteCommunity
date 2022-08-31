#ifndef LOUVAIN_H
#define LOUVAIN_H

#include "../optimiser.h"

using namespace std; 

template <typename T>
class Louvain : public Optimiser<T> {
	public:
	
		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;

		//functions
		using Optimiser<T>::greedy_update; 
		using Optimiser<T>::assign_labels; 
		using Optimiser<T>::setup_sizes; 


		Louvain(T& Q);	
		double optimise();

};

template <typename T>
Louvain<T>::Louvain(T& _Q) : Optimiser<T>(_Q) {
		optimiser_type = "louvain";
}

template <typename T>
double Louvain<T>::optimise() {
	
	double dq = 0;
	int level = 0;
	vector< unordered_map< int, int > > community_maps;

	do{

		//find new community labels	
		dq = greedy_update(); 

		//save community assignment
		community_maps.push_back( Q.node_to_comm );
		
		//use communities to make induced graph
		WeightedNetwork GI = Q.G.induced_graph( Q.node_to_comm ); 	

		
		//reset quality function with GI
		Q.setup_network(GI);
		setup_sizes();
		
		++level;
		
	} while(dq > eps);

	Q.eval();	
	assign_labels(community_maps);	
	
	return 0;
}


#endif
