#ifndef OPTIMISER_H
#define OPTIMISER_H

#include "../qualities/quality.h"
#include "../qualities/modularity.h"
#include "../qualities/barber_modularity.h"
#include "../qualities/projected_modularity.h"
#include "../qualities/cp_modularity.h"
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>

using namespace std; 

//interface for community optimiser
template <typename T>
class Optimiser {
	public:
		
		T &Q; 
		string optimiser_type;
		double val;
		int num_nodes;
		int num_communities;
		double eps;
		unordered_map<int, int> labels;
		int passes;
	
		void setup_sizes();
		Optimiser(T& _Q);

		virtual double optimise() = 0;
		void print_basic();
		void shuffle(vector<int> &x);		//randomise the node order
		double greedy_update(); 			//do one step of louvain
		void fixupQ();

		void print_labels();
		void assign_labels(vector< unordered_map< int, int > > community_maps); //unfold the louvain algo
};


template <typename T> Optimiser<T>::Optimiser(T& _Q) : Q(_Q) {
	
	optimiser_type = "optimiser";
	for(auto &i : Q.G.get_nbrs){ labels[ i.first ] = i.first; }
	setup_sizes();
	eps = 1e-8;
	
}
template <typename T> void Optimiser<T>::setup_sizes(){
	num_nodes = Q.G.num_nodes;
	num_communities = Q.num_communities;
	
}

template <typename T> void Optimiser<T>::print_basic(){
	cout << "This is a " << optimiser_type << " optimiser" << endl;
	cout << "Completed " << passes << " passes" << endl;
	cout << "Last computed " << Q.quality_type << " = " << Q.val << endl;	
	print_labels();
}

///shuffle the elements of x in place
template <typename T> void Optimiser<T>::shuffle(vector<int> &x){
	
	int n = x.size();
	for(int i=0; i<n; ++i){
		int ridx = rand()%n;
		int tmp = x[i];
		x[i] = x[ridx];
		x[ridx] = tmp;
	}
	
}

///Equivalent to one_level of a Louvain typr algorithm
template <typename T> double Optimiser<T>::greedy_update(){
	
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

			//random neighbour visit order
			vector<int> nbr_list; for(auto &i : Q.G.get_nbrs[node]){ nbr_list.push_back(i.first); } shuffle(nbr_list);
						
			for(auto &nbr : nbr_list){ 
				if( node == nbr ){ continue; } //self or already in same comunity as nbr, dq=0
				
				double g = Q.gain(node, Q.node_to_comm[nbr]); //calculate quality change from adding 'node' to community of 'nbr'

				if( g > dq ){ //best change seen so far
					best_comm = Q.node_to_comm[nbr];
					dq = g;
				} //should put in a tie breaker
			}
			
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

///Q fix inconsistencies in Q caused by node swapping
template <typename T> void Optimiser<T>::fixupQ(){
	Q.setup_comm_to_node();
	Q.init();
	Q.eval();
}



///compute final communtiy assignment after louvain steps
template <typename T> void Optimiser<T>::assign_labels(vector< unordered_map< int, int > > community_maps){

	//community assignment
	for(unsigned i=0; i<community_maps.size(); ++i){
		for(auto &j : labels){
			labels[ j.first ] = community_maps[i][ j.second ];
		}
	}
	
}

template <typename T> void Optimiser<T>::print_labels(){
	for(auto &i : labels){ cout << "Label (" << i.first << ") : " << i.second << endl; }		
}











#endif
