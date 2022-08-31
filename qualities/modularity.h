#ifndef MODULARITY_H
#define MODULARITY_H

#include "quality.h"

using namespace std; 

/**
 * Classic Newman modularity for the standard Louvain algorithm
 * 
 * 
 * Q = 1/2E sum_ij (A_ij - ki kj/2E) delta( c(i) , c(j) ) 
 * 
 * sum_ij A_ij delta( c(i) , c(j) ) = in[ c ] = sum of edges in community
 * sum_ij (ki kj/2E) delta( c(i) , c(j) ) = tot[ c ]^2. tot[ c ] = sum of degrees in community
 * 
 * Q = (1/2E) * sum_c in[c] - tot[c]^2/2E 
 * 
 * remove/add node n from c:
 * 	in[c] -/+ sum (links from n to other nodes in c) + Ann = 2knc + Ann
 *  tot[c] -/+ degree(n) = kn
 * 
 * add node n from isolated communtiy to c
 * dQ = (1/2E) * (  ((in[c]+2knc+Ann) - (in[c]+Ann)) -
 * 					((tot[c] + kn)**2 - (tot[c]**2 + kn**2))/2E  )
 * 
 *    = (1/2E) * ( 2knc - 2*tot[c]*kn )
 * 
 * merge c1 and c2
 * in[c1] = in[c1]+in[c2] + (edges between c1 and c2) = in[c1]+in[c2] + 2*v12
 * tot[c1] = tot[c1] + tot[c2]
 * 
 * dQ = (1/2E) * ( (in[c1]+in[c2] + 2*v12 - in[c1]+in[c2]) -
 * 				   ( (tot[c1]**2 + 2tot[c1]tot[c2] + tot[c2]**2) - (tot[c1]**2 + tot[c2]**2))/2E  ) 
 *    = (1/2E) * ( 2*v12 - 2*tot[c1]tot[c2]/2E )
 *
 * 
 * swap node from c1 to c2.
 * 
 * dQ = (1/2E) * (  (in[c2] + 2knc2 + Ann - in[c1] - 2knc2 - Ann) -
 * 				    ( (tot[c2]+kn)^2 + (tot[c1]-kn)^2 - tot[c1]^2 - tot[c2]^2)/2E ) 
 *    = (1/2E) * ( 2(knc2 - knc1) - 2kn(tot[c2] - tot[c1] + kn)/2E )
 **/
template<typename T>
class Modularity : public Quality<T> {
	public:
	
		//This is necessary for a templated class
		//Doing this way avoids having to use this-> in the method definitions
		using Quality<T>::num_communities;
		using Quality<T>::val;
		using Quality<T>::gamma;
		using Quality<T>::min;
		using Quality<T>::norm;
		using Quality<T>::node_to_comm;
		using Quality<T>::comm_to_node;
		using Quality<T>::in;
		using Quality<T>::tot;
		using Quality<T>::G; 
		using Quality<T>::quality_type;
		using Quality<T>::dc;
		using Quality<T>::v12;
		using Quality<T>::swap;
		using Quality<T>::update_nodemaps;

		Modularity(T& _G);
		Modularity(T& _G, unordered_map<int, int> &labels);
	
		void init();
		void eval();
		void set_norm();
		void remove(int node, int comm);
		void insert(int node, int comm);
		void join(int comm1, int comm2);

		double gain(int node, int comm);
		double swap_gain(int node, int comm1, int comm2);
		double join_gain(int comm1, int comm2);
};



template <typename T>
Modularity<T>::Modularity(T& _G) : Quality<T>(_G){
	quality_type = "modularity";	
	Modularity<T>::init();
}

template <typename T>
Modularity<T>::Modularity(T& _G, unordered_map<int, int> &labels) : Quality<T>(_G, labels){
	quality_type = "modularity";	
	Modularity<T>::init();
}

template <typename T>
void Modularity<T>::set_norm(){
	norm = 2/G.num_links;
}

///set up in and tot vectors
template <typename T>
void Modularity<T>::init(){

	in.clear();
	tot.clear();
	for(auto &i : comm_to_node){ in[i.first] = 0; tot[i.first] = 0; }
	
	//loop over edges
	for(auto &i: G.get_nbrs){
		for(auto &j: i.second){
			if( node_to_comm[i.first] == node_to_comm[ j.first ] ){
				in[ node_to_comm[i.first] ] += j.second;
			} 
		}
	}

	//loop over nodes
	for(auto &i : G.degrees){ tot[ node_to_comm[i.first] ] += i.second; }
}

template <typename T>
void Modularity<T>::eval(){
	
	double q = 0;
	for(auto &i : comm_to_node){
		if(i.second.size()>0){		//skip empty communities
			q += in[i.first] - gamma*tot[i.first]*tot[i.first]/G.num_links;
		}
	}
	val = min*q/G.num_links;
}

template <typename T>
void Modularity<T>::remove(int node, int comm){

	in[comm]  -= 2*dc[comm] + G.selfloop(node);
	tot[comm] -= G.degrees[node];

	update_nodemaps(node, comm, -1);
}

template <typename T>
void Modularity<T>::insert(int node, int comm){

	in[comm]  += 2*dc[comm] + G.selfloop(node);
	tot[comm] += G.degrees[node];

	update_nodemaps(node, -1, comm);
}



template <typename T>
double Modularity<T>::gain(int node, int comm) {

	return min*norm*(dc[comm] - gamma*tot[comm]*G.degrees[node]/G.num_links); 
}

//node goes from comm1 to comm2
template <typename T>
double Modularity<T>::swap_gain(int node, int comm1, int comm2) {

	return min*norm*((dc[comm2] - dc[comm1]) - gamma*(tot[comm2] - tot[comm1] + G.degrees[node])*G.degrees[node]/G.num_links); 
}

template <typename T>
void Modularity<T>::join(int comm1, int comm2){


	in[comm2]  += in[comm1] + 2*v12;
	in[comm1] = 0;

	tot[comm2] += tot[comm1];
	tot[comm1] = 0;
	
	Modularity<T>::qjoin(comm1, comm2);
  
}


template <typename T>
double Modularity<T>::join_gain(int comm1, int comm2) {

	return min*norm*(v12 - gamma*tot[comm1]*tot[comm2]/G.num_links); 

}



#endif
