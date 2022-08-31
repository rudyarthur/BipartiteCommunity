#ifndef PROJECTEDLOUVAIN_H
#define PROJECTEDLOUVAIN_H


#include "../optimiser.h"

using namespace std; 

template <typename T>
class ProjectedLouvain : public Optimiser< T > {
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
					
		ProjectedLouvain(T& _Q);		
		double optimise();
};

template <typename T>
ProjectedLouvain<T>::ProjectedLouvain(T& _Q) : Optimiser<T>(_Q) {
		if(Q.quality_type != "projected_modularity"){
			cerr << "only use this optimiser with projected modularity!" << endl;
			exit(1);
		}
		optimiser_type = "projected_louvain";
}		


template <typename T>
double ProjectedLouvain<T>::optimise() {

	double dq = 0;
	int level = 0;
	vector< unordered_map< int, int > > community_maps;

	do{

		//find new community labels	
		dq = greedy_update(); 
		
		
		//save community assignment
		community_maps.push_back( Q.node_to_comm );
		
		//use communities to make induced graph
		//With include loops=0 this is unnecessary otherwise we lose edges weight as we aggregate
		//GP is the projection (unipartite). GI uses comm labels to induce a smaller graph
		WeightedNetwork GI = Q.P.GP.induced_graph( Q.node_to_comm ); 
		//G is the bipartite network. BI induces the graph on the correct size of the bipartite
		BipartiteNetwork BI = Q.P.G.induced_graph( Q.node_to_comm, Q.P.side);
		

		//Reset Quality	
		Q.P.GP = GI; //GI is NOT the projection of BI in general
		Q.P.G = BI; 
		//make the new Q
		Q.setup_network(GI);
		setup_sizes();

		++level;
		
	} while(dq > eps);

	Q.eval();	
	assign_labels(community_maps);	
	
	return 0;
}


#endif
