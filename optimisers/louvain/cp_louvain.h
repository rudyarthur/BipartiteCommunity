#ifndef CPLOUVAIN_H
#define CPLOUVAIN_H

#include "../optimiser.h"

using namespace std; 

//template <typename T>
class CPLouvain : public Optimiser< CPModularity<WeightedNetwork> > {
	public:

		CPLouvain( CPModularity<WeightedNetwork>& Q);	
		double optimise();
		double greedy_update();
		
};

//template <typename T>
CPLouvain::CPLouvain(CPModularity<WeightedNetwork>& _Q) : Optimiser(_Q) {
		optimiser_type = "cp_louvain";
}

///Equivalent to one_level of a Louvain typr algorithm
double CPLouvain::greedy_update(){
	
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
			
			int ci = comm/2;
			int xi = comm%2; 	
			for(int xi = 0; xi<2; ++xi){ //same community, different parity
				double g = Q.gain(node,2*ci+(xi)); //calculate quality change from adding 'node' to community of 'nbr'

				if( g > dq ){ //best change seen so far
					best_comm = 2*ci+(xi);
					dq = g;
				}
			}
			for(auto &nbr : nbr_list){ 
				if( node == nbr ){ continue; } //self
				int cj = Q.node_to_comm[nbr]/2;
				for(int xj=0; xj<2; ++xj){ //both parities
					double g = Q.gain(node, 2*cj + xj); //calculate quality change from adding 'node' to community of 'nbr'
					if( g > dq ){ //best change seen so far
						best_comm = 2*cj + xj;
						dq = g;
					} 
				}
			}
			//push node to comm
 			Q.insert(node, best_comm);

 			//check if something changed
			if (best_comm!=comm){++moves;}
			
		} //finish node loop

  //if we moved a node we must have improved the quality
  } while (moves>0);
  
  fixupQ();
  
  return Q.val-old_val;
 
	
}

double CPLouvain::optimise() {
	
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
		Q.node_to_cp.clear(); for(auto &i : Q.comm_to_node){ Q.node_to_cp[i.first/2] = i.first%2; }
		Q.setup_network(GI);
		setup_sizes();
		
		++level;
		
	} while(dq > eps);

	Q.eval();	
	assign_labels(community_maps);	
	
	return 0;

}


#endif
