#ifndef INVLOUVAIN_H
#define INVLOUVAIN_H

#include "../optimiser.h"
#include <iomanip>

using namespace std; 

/**
 * An interesting method for finding hidden bipartite structure 
 * https://www.sciencedirect.com/science/article/pii/S1389128621000621
 * **/
template <typename T>
class InvLouvain : public Optimiser<T> {
	public:
			//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::passes;

		//functions
		using Optimiser<T>::greedy_update; 
		using Optimiser<T>::assign_labels; 
		using Optimiser<T>::setup_sizes;
		using Optimiser<T>::shuffle;
		using Optimiser<T>::fixupQ;
		
		
		InvLouvain(T& Q);	
		double optimise();
		double greedy_update();
};

template <typename T>
InvLouvain<T>::InvLouvain(T& _Q) : Optimiser<T>(_Q) {
		optimiser_type = "inv_louvain";
		Q.set_min();
}

//It's a bit annoying to have to re-implement all of this just because we're doing next nearest neighbours...
template <typename T> double InvLouvain<T>::greedy_update(){
	
	double dq;
	int moves;
	passes = 0;
    double old_val = Q.val;
	
		
	do {
		moves = 0;
		//random node visit order
		vector<int> node_list; for(auto &i : Q.node_to_comm){ node_list.push_back(i.first); } shuffle(node_list); //random order
		
		for (auto &node : node_list) {

			int comm = Q.node_to_comm[node];
						

			//get links from node to other communities
			Q.links_to_comms(node);
			//pop node from comm
			Q.remove(node, comm);
			int best_comm = comm;
			dq = 0;		


			//random next neighbour visit order
			unordered_set<int> nbr_set; 
			for(auto &i : Q.G.get_nbrs[node]){
				for(auto &j : Q.G.get_nbrs[ i.first ]){ 
					nbr_set.insert( j.first ); 
				}
			} 
			vector<int> nbr_list; for(auto &i : nbr_set){nbr_list.push_back(i);} shuffle(nbr_list);
			

			for(auto &nbr : nbr_list){ 
				if( node == nbr ){ continue; } //self or already in same comunity as nbr, dq=0
				
				double g = Q.gain(node, Q.node_to_comm[nbr]); //calculate quality change from adding 'node' to community of 'nbr'
				if( g > dq ){ //best change seen so far, overwrite the others
					best_comm = Q.node_to_comm[nbr];
					dq = g;
				} //TODO should put in a tie breaker
			
			}//TODO what if keeping node in separate community is optimal?
			
			//push node to comm
 			Q.insert(node, best_comm);
 			//check if something changed
			if (best_comm!=comm){++moves;}
			
		} //finish node loop
		++passes;
		
  //if we moved a node we must have improved the quality
  } while (moves>0);
  
  fixupQ();
  
  return Q.val-old_val;
	
}


template <typename T>
double InvLouvain<T>::optimise() {
	
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
