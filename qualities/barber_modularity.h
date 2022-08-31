#ifndef BARBERMODULARITY_H
#define BARBERMODULARITY_H

#include "quality.h"

using namespace std; 

/**
 * Barber modularity for the standard Louvain algorithm
 * 
 * F = sum_ij B_ij
 * qi = sum_k Bik
 * dk = sum_i Bik
 * 
 * 
 * 
 * Q = 1/F sum_ij (B_ij - qi dj /F) delta( c(i) , c(j) ) 
 * 
 * sum_ij B_ij delta( c(i) , c(j) ) = in[ c ] = sum of edges in community
 * sum_ij (qi dj /F) delta( c(i) , c(j) ) = totq[ c ] totd[ c ] = (sum of left degrees in community)  x (sum of right degrees in community)
 * 
 * Q = (1/F) * sum_c (in[c] - totq[c]totd[c]/F) 
 * 
 * remove/add left node n from c:
 * 	in[c] -/+ sum (links from n to other nodes in c) = qnc
 *  totq[c] -/+ qn
 * 
 * left->right  q->d
 * 
 * add node n from isolated communtiy to c
 * dQ = (1/F) * (  ((in[c]+qnc) - (in[c])) -
 * 					((totq[c] + qn)totd[c] - (totq[c]totd[c]))/F  )
 * 
 *    = (1/F) * ( qnc - qn totd[c]/F )
 * 
 * merge c1 and c2
 * in[c1] = in[c1]+in[c2] + (edges between c1 and c2) = in[c1]+in[c2] + v12
 * totq[c1] = totq[c1] + totq[c2]
 * 
 * dQ = (1/F) * ( (in[c1]+in[c2] + v12 - in[c1]-in[c2]) -
 * 				   ( (totq[c1] + totq[c2])(totd[c1] + totd[c2]) - (totq[c1]totd[c1] + totq[c2]totd[c2]))/F  ) 
 *    = (1/F) * ( v12 - (totq[c1]totd[c2] + totq[c2]*totd[c1])/F )
 *
 * 
 * swap node from c1 to c2.
 * 
 * dQ = (1/F) * (  (in[c2] + qnc2 + in[c1] - qnc1) - (in[c2] + in[c1])-
 * 				    ( (totq[c2]+qn)totd[c2] + (totq[c1]-qn)totd[c1] - totq[c2]totd[c2] - totq[c1]totd[c1])/F ) 
 *    = (1/F) * ( (qnc2 - qnc1) - qn(totd[c2] - totd[c1])/F )
 **/
class BarberModularity : public Quality< BipartiteNetwork > {
	public:
		
		BarberModularity(BipartiteNetwork& _G);
		BarberModularity(BipartiteNetwork& _G, unordered_map<int, int> &labels);

		unordered_map<int, double> totr;	//the degrees on the other side
		double F;
		
		void init();
		void eval();
		void set_norm();		
		void remove(int node, int comm);
		void insert(int node, int comm);
		void join(int comm1, int comm2);

		double gain(int node, int comm);
		double swap_gain(int node, int comm1, int comm2);
		double join_gain(int comm1, int comm2);
		
		void print_basic();
		void print_partial_sums();
};


BarberModularity::BarberModularity(BipartiteNetwork& _G) : Quality<BipartiteNetwork>(_G){
	quality_type = "barber_modularity";	
	BarberModularity::init();
}

BarberModularity::BarberModularity(BipartiteNetwork& _G, unordered_map<int, int> &labels) : Quality<BipartiteNetwork>(_G, labels){
	quality_type = "barber_modularity";	
	BarberModularity::init();
}

void BarberModularity::set_norm(){
	norm = 1.0/F;
}

///output basic info
void BarberModularity::print_basic(){
	Quality<BipartiteNetwork>::print_basic();
	cout << "F = " << F << endl;
}

///print partial sums
void BarberModularity::print_partial_sums(){
	for(auto &i : in){ cout << " in[" << i.first << "] = " << i.second << " tot[" << i.first << "] = " << tot[i.first] << " totr[" << i.first << "] = " << totr[i.first] << endl; }
}

void BarberModularity::init(){

	in.clear();
	tot.clear();
	totr.clear();
	for(auto &i : comm_to_node){ in[i.first] = 0; tot[i.first] = 0; totr[i.first] = 0; }
	
	//loop over edges
	for(auto &i: G.get_nbrs){
		for(auto &j: i.second){
			if( node_to_comm[i.first] == node_to_comm[ j.first ] ){
				in[ node_to_comm[i.first] ] += j.second*0.5; //not double counting here
			} 
		}
	}

	//loop over nodes
	for(auto &i : G.degrees){ 
		if(G.colour[i.first] == 0){
			tot[ node_to_comm[i.first] ] += i.second; //left degree sums
		} else {
			totr[ node_to_comm[i.first] ] += i.second; //right degree sums
		}
	}
	F = G.num_links/2.0;	  //factor of 2 is for double counting edges
	
}


void BarberModularity::eval(){
	
	double q = 0;
	for(auto &i : comm_to_node){
		if(i.second.size()>0){		//skip empty communities
			q += in[i.first] - gamma*tot[i.first]*totr[i.first]/F;
		}
	}
	val = min*q/F;

}


void BarberModularity::remove(int node, int comm){


	in[comm]  -= dc[comm];
	if( G.colour[node] == 0 ){
		tot[comm] -= G.degrees[node];
	} else {
		totr[comm] -= G.degrees[node];
	}
	
	update_nodemaps(node, comm, -1);

}


void BarberModularity::insert(int node, int comm){

	in[comm]  += dc[comm];
	if( G.colour[node] == 0 ){
		tot[comm] += G.degrees[node];
	} else {
		totr[comm] += G.degrees[node];
	}
	
	update_nodemaps(node, -1, comm);

}

double BarberModularity::gain(int node, int comm) {


	if( G.colour[node] == 0 ){
		return min*norm*(dc[comm] - gamma*totr[comm]*G.degrees[node]/F);
	} 
	return min*norm*(dc[comm] - gamma*tot[comm]*G.degrees[node]/F);

	
}

double BarberModularity::swap_gain(int node, int comm1, int comm2) {
	if( G.colour[node] == 0 ){
		return min*norm*((dc[comm2] - dc[comm1]) - gamma*(totr[comm2] - totr[comm1])*G.degrees[node]/F); 
	}
	return min*norm*((dc[comm2] - dc[comm1]) - gamma*(tot[comm2] - tot[comm1])*G.degrees[node]/F); 

}

void BarberModularity::join(int comm1, int comm2){

	in[comm2]  += in[comm1] + v12;
	in[comm1] = 0;

	tot[comm2] += tot[comm1];
	tot[comm1] = 0;

	totr[comm2] += totr[comm1];
	totr[comm1] = 0;
		
	BarberModularity::qjoin(comm1, comm2);
	
}

double BarberModularity::join_gain(int comm1, int comm2) {
	 
	return min*norm*(v12 - gamma*(tot[comm1]*totr[comm2] + tot[comm2]*totr[comm1])/F); 

}



#endif
