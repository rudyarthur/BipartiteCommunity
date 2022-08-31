#ifndef PROJECTEDMODULARITY_H
#define PROJECTEDMODULARITY_H

#include "quality.h"

using namespace std; 

/**
 * Projected modularity for the standard Louvain algorithm
 * 
 * F = sum_ij B_ij
 * qi = sum_k Bik
 * dk = sum_i Bik
 * A = BB^T
 * 2E = sum_ij A_ij
 * ki = sum_j A_ij
 * 
 * 
 * Q = 1/2E sum_ij (A_ij - D2(qi qj /F(F-s)) ) delta( c(i) , c(j) ) 
 * N = D2/F(F-s)
 * 
 * sum_ij A_ij delta( c(i) , c(j) ) = in[ c ] = sum of edges in community
 * sum_ij qi qj delta( c(i) , c(j) ) = totq[ c ]^2 = (sum of left degrees in community)^2
 * 
 * Q = (1/2E) * sum_c (in[c] - N*totq[c]^2 )
 * 
 * remove/add left node n from c:
 * 	in[c] -/+ sum (links from n to other nodes in c) = 2knc + Ann
 *  totq[c] -/+ qn
 * 
 * left->right  q->d
 * 
 * add node n from isolated communtiy to c
 * dQ = (1/2E) * (  ((in[c]+2knc + Ann) - (in[c]+Ann)  
 * 					((totq[c] + qn)^2 - (totq[c]^2 + qn^2))N  )
 * 
 *    = (1/2E) * ( 2knc - (2qn totq[c])N )
 * 
 * merge c1 and c2
 * in[c1] = in[c1]+in[c2] + (edges between c1 and c2) = in[c1]+in[c2] + 2v12
 * totq[c1] = totq[c1] + totq[c2]
 * 
 * dQ = (1/2E) * ( (in[c1]+in[c2] + 2*v12 - in[c1]+in[c2]) -
 * 				   ( (totq[c1]**2 + 2totq[c1]totq[c2] + totq[c2]**2) - (totq[c1]**2 + totq[c2]**2))*N  ) 
 *    = (1/2E) * ( 2*v12 - 2*totq[c1]totq[c2]*N )
 *
 * 
 * swap node from c1 to c2.
 * 
 * dQ = (1/2E) * (  (in[c2] + 2knc2 + Ann - in[c1] - 2knc2 - Ann) -
 * 				    ( (totq[c2]+qn)^2 + (totq[c1]-qn)^2 - totq[c1]^2 - totq[c2]^2)*N ) 
 *    = (1/2E) * ( 2(knc2 - knc1) - 2qn(totq[c2] - totq[c1] + qn)*N )
 **/
template <typename T>
class ProjectedModularity : public Quality< WeightedNetwork > {
	public:
	
		T &P;		//this is going to be a projector
		double F;	//Edges in bipartite
		double N;	//The normalization factor
		
		ProjectedModularity(T &_p);
		ProjectedModularity(T &_p, WeightedNetwork& _G);
		ProjectedModularity(T &_p, WeightedNetwork& _G, unordered_map<int, int> &labels);
		ProjectedModularity(T &_p, unordered_map<int, int> &labels);
		
			
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
ProjectedModularity<T>::ProjectedModularity(T &_p) :  Quality< WeightedNetwork >(_p.GP), P(_p)  {
	quality_type = "projected_modularity";	
	ProjectedModularity<T>::init();
}
template <typename T>
ProjectedModularity<T>::ProjectedModularity(T &_p, unordered_map<int, int> &labels) :  Quality< WeightedNetwork >(_p.GP, labels), P(_p)  {
	quality_type = "projected_modularity";	
	ProjectedModularity<T>::init();
}
template <typename T>
ProjectedModularity<T>::ProjectedModularity(T &_p, WeightedNetwork& _G) :  Quality< WeightedNetwork >(_G), P(_p)  {
	quality_type = "projected_modularity";	
	ProjectedModularity<T>::init();
}

template <typename T>
ProjectedModularity<T>::ProjectedModularity(T &_p, WeightedNetwork& _G, unordered_map<int, int> &labels) : Quality< WeightedNetwork >(_G, labels), P(_p)  {
	quality_type = "projected_modularity";	
	if(labels.size() != (unsigned)G.num_nodes){ cerr << "need one label per node in projected network for " << quality_type << endl; exit(1); }
	ProjectedModularity<T>::init();
}

template <typename T>
void ProjectedModularity<T>::set_norm(){
	norm = 2.0/G.num_links;
}

template <typename T>
void ProjectedModularity<T>::init(){

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
	for(auto &i : P.G.degrees){ 
		if(P.G.colour[i.first] == P.side){
			tot[ node_to_comm[i.first] ] += i.second; //bipartite degree sums
		}
	}
	F = P.G.num_links/2.0;	  //factor of 2 is for double counting edges
	N = P.D2/(F*(F-P.inc));

}


template <typename T>
void ProjectedModularity<T>::eval(){
	
	double q = 0;
	for(auto &i : comm_to_node){
		if(i.second.size()>0){		//skip empty communities
			q += in[i.first] - gamma*tot[i.first]*tot[i.first]*N;	
			//cout << "q += " << 		in[i.first] << " - " << tot[i.first] << "**2" << endl;
		}
	}
	val = min*q/G.num_links;
	

}

template <typename T>
void ProjectedModularity<T>::remove(int node, int comm){

	in[comm]  -= 2*dc[comm] + G.selfloop(node);
	tot[comm] -= P.G.degrees[node];

	update_nodemaps(node, comm, -1);

}

template <typename T>
void ProjectedModularity<T>::insert(int node, int comm){

	in[comm]  += 2*dc[comm] + G.selfloop(node);
	tot[comm] += P.G.degrees[node];

	update_nodemaps(node, -1, comm);

}

template <typename T>
double ProjectedModularity<T>::gain(int node, int comm) {
	
	return min*norm*(dc[comm] - gamma*tot[comm]*P.G.degrees[node]*N); 
	
}

template <typename T>
double ProjectedModularity<T>::swap_gain(int node, int comm1, int comm2) {
	
	return min*norm*((dc[comm2] - dc[comm1]) - gamma*(tot[comm2] - tot[comm1] + P.G.degrees[node])*P.G.degrees[node]*N  );

}

template <typename T>
void ProjectedModularity<T>::join(int comm1, int comm2){


	in[comm2]  += in[comm1] + 2*v12;
	in[comm1] = 0;

	tot[comm2] += tot[comm1];
	tot[comm1] = 0;
	
	ProjectedModularity<T>::qjoin(comm1, comm2);
  
}



template <typename T>
double ProjectedModularity<T>::join_gain(int comm1, int comm2) {

	return min*norm*(v12 - gamma*(tot[comm1]*tot[comm2])*N );

}

#endif
