#ifndef CPMODULARITY_H
#define CPMODULARITY_H

#include "quality.h"

using namespace std; 

/**
 * CP modularity for the standard Louvain algorithm
 * 
 * xi = CP label
 * Q = 1/2E sum_ij (A_ij - ki kj/2E) (xi + xj - xi*xj) delta( c(i) , c(j) ) 
 * 
 * sum_ij A_ij (xi + xj - xi*xj) delta( c(i) , c(j) ) = in[ c ] 
 * sum_ij (ki kj/2E) (xi + xj - xi*xj) delta( c(i) , c(j) ) = totc[ c ]*( 2*tot[c] - totc[c] ) 
 * 
 * 
 * 
 * in[c] = sum_{ij in c} A_ij (xi + xj - xi*xj) 
 * in[c+n] =  sum_{ij in c+n} A_ij (xi + xj - xi*xj) 
 *         =  sum_{ij in c} A_ij (xi + xj - xi*xj) + 2*sum_{i in c} A_in (xi + xn - xi*xn) +  A_nn (2xn - xn**2) //0**2 = 0, 1**2 = 1
 *         =  sum_{ij in c} A_ij (xi + xj - xi*xj) + 2*sum_{i in c} A_in (xi + xn - xi*xn) +  A_nn (xn)
 *         =  in[c] + 2*sum_{i in c} A_in (xi + xn - xi*xn) +  A_nn (xn)
 * 
 * inp[c] = sum_{ij in c} A_ij xi
 * inp[c+n] = inp[c] + sum_{i in c} A_in xn + A_nn xn = inp + dc[n] + Ann xn
 * dc[n] = sum_{i in c} A_in xn = xn sum_{i in c} A_in = xn kc
 * 
 * inc[c] = sum_{ij in c} A_ij xi xj
 * inc[c+n] = inc[c] + sum_{i in c} A_in xi xn + Ann xn = inc + dcc[n] + Ann xn
 * dcc[c] = sum_{i in c} A_in xi xn = xn sum_{i in c} A_in xi
 * 
 * sum_ij (ki kj) (xi + xj - xi*xj) delta( c(i) , c(j) ) = 2 sum_{ij in c} ki kj xi - sum_{ij in c} ki xi kj xj
 * 														 = 2 tot * totc - totc^2
 * 														 = totc( 2 tot - totc )
 * tot[c] = sum_{i in c} k_i
 * tot[c+n] = tot[c] + kn
 * totc[a] = sum_{i in a} ki xi
 * totc[a+n] = totc[a] + xn kn
 * 
 * Q = (1/2E) * sum_c ( 2*inp [c] - inc [c]  - totc[ c ]*( 2*tot[c] - totc[c] ) /2E  )
 * 
 * 
 * add node n from isolated communtiy to c as xn
 *in_diff = 2*(inp + xndc[n] + Ann xn) - (inc + xndcc[n] + Ann xn)
 * 		  - (2*inp - inc + Ann)
 * 		  = 2 xdc + 2Ann x - xdcc - Ann x - Ann
 * 		  = xn(2dc[n] - dcc[n]) - Ann (1-xn)
 * 					
 * tot_diff = (totc + xk)*(2tot + 2k - totc - xk)  - totc*( 2*tot - totc ) - k^2
 *          = totc*(2tot - totc)  + totc(2k - xk) + xk*(2*tot + 2k - totc - xk) - totc*( 2*tot - totc ) - k^2
 * 			= totc(2k - xk - xk) + tot(2xk) + k^2(2x - x^2 - 1)
 * 			= 2totc(k(1-x)) + tot(2xk) - k^2(1-x)^2
 * 			= (1-x)*( k*(2totc - k) ) + (x) * (tot * 2k)
 * 
 * swap node from c1 to c2. where it is xn1 in c1 and xn2 in c2
 * 
 * in part :   2*(inp[c2] + x2dc2[n] + Ann xn2) - (inc[c2] + x2dcc2[n] + Ann xn2) + 2*(inp[c1] - x1dc1[n] - Ann xn1) - (inc[c1] - x1dcc1[n] - Ann xn1)
 * 			 -(2*x2(inp[c2] - inc[c2]) + 2*x1(inp[c1] - inc[c1])) 
 * 			= 2(x2dc2 - x1dc1) - (x2dcc2 - x1dcc1) + Ann(x2 - x1)
 * 
 * totc[a+n] = totc[a] + xn kn
 * tot[c+n] = tot[c] + kn
 * 
 * 
 * tot part : (totc[c2] + x2 k)*( 2tot[c2] + 2k - totc[c2] - x2 k ) +  (totc[c1] - x1 k)*( 2tot[c1] - 2k - totc[c1] + x1 k )
 *           -( totc[c2]*( 2tot[c2] - totc[c2] ) + totc[c1]*( 2tot[c1] - totc[c1] ) )
 * 			= totc[c2](2k-x2 k) + x2 k(2tot[c2] + 2k - totc[c2] - x2 k) - totc[c1](2k - x1 k) - x1 k( 2tot[c1] - 2k - totc[c1] + x1 k )
 * 			=  + x2 k( + 2k  - x2 k) - x1 k(  - 2k + x1 k )
 *          = 2k( totc2 (1-x2) - totc1 (1-x1) ) + 2k (tot2 x2 - tot1 x1) + k^2( x2(2-x2) + x1(2-x1) )
 *          = k*( 2( totc2 (1-x2) - totc1 (1-x1) ) + 2 (tot2 x2 - tot1 x1) + k( x2(2-x2) + x1(2-x1) ) )
 * 
 * 
 * 	
 * 
 * merge c1x1 and c2x2 (e.g. core of c1 goes into the periphery of c2)
 * inp[c1] = sum_{ij in c1} A_ij xi
 * inp[c2] = sum_{ij in c2} A_ij xi													     
 * inp[c1+c2] = inp[c2] + (join new nodes to c2 core)	     + (if adding as core)   (add core-core)
 * inp[c1+c2] = inp[c2] + sum_{i in c1x1, j in c2} A_ij xj + x2                  * (sum_{i in c1x1, j in c1x1} A_ij)
 * v12 = sum_{i in c1x1, j in c2} A_ij xj 
 * vp1 = sum_{i in c1x1, j in c1x1} A_ij = x1*in[c1] + (don't accumuate the pp sum, so have to calculate it)
 * inp[c1+c2] = inp[c2] + vp12 + x2*vp1
 * inp[c1] = (1-x1)*(inc[c1])
 * 
 * inc[c1] = sum_{ij in c1} A_ij xi xj
 * inc[c2] = sum_{ij in c2} A_ij xi xj						
 * inc[c1+c2] = inc[c2] + (if core)         (new nodes to old core)		+  (new core to itself)	
 * inc[c1+c2] = inc[c2] + x2*       (sum_{i in c1x1, j in c2} A_ij xj + (sum_{i in c1x1, j in c1x1} A_ij))
 * v12 = sum_{i in c1x1, j in c2} A_ij xj
 * inc[c1+c2] = inc[c2] + x2*( v12 + v1 )
 * inc[c1] = (1-x1)*(inc[c1])
 * 
 * 
 * 
 * tot[c1] = sum_{i in c1} k_i
 * tot[c2] = sum_{i in c2} k_i
 * tot[c1x1+c2] = tot[c2] + x1*totc[c1] + (1-x1)*(tot[c1] - totc[c1]) 
 * 				= tot[c2] + (tot[c1]-totc[c1]) + x1*(2totc[c1] - tot[c1])
 * tot[c1] = tot[c1] - x1*totc[c1] - (1-x1)*(tot[c1] - totc[c1])
 * 		   = tot[c1] -(tot[c1] - totc[c1]) - x1(2totc[c1] - tot[c1])
 * 
 * totc[c1] = sum_{i in c1} ki xi
 * totc[c2] = sum_{i in c2} ki xi
 * totc[c1x1+c2] = totc[c2] + x2*(x1*totc[c1] + (1-x1)*(tot[c1] - totc[c1]))
 *               = totc[c2] + x2*(  (tot[c1] - totc[c1]) + x1*(2totc[c1] - tot[c1])  )
 * totc[c1] = totc[c1] - x1*totc[c1]
 * 
 * 
 * dQ = too complicated...
 * 			
 * 
 **/
template<typename T>
class CPModularity : public Quality<T> {
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
		using Quality<T>::vc1;
		using Quality<T>::swap;
		using Quality<T>::update_nodemaps;
		
		//1 is core 0 is periphery
		unordered_map<int, int> node_to_cp;
		unordered_map<int, double> totc;	//sum of core degrees
		unordered_map<int, double> inc;		//sum of core edges
		unordered_map<int, double> dcc; 	//links to neighbour communities
		
		
		CPModularity(T& _G);
		CPModularity(T& _G, unordered_map<int, int> &labels);
		CPModularity(T& _G, unordered_map<int, int> &labels, unordered_map<int, int> &cp_labels);

		void init();
		void eval();
		void set_norm();
		void links_to_comms(int node);	//calculate the connections from node to all the comms		
		void links_between(int comm1, int comm2); //calculate the number of links between c1 and c2
		
		void remove(int node, int comm);
		void insert(int node, int comm);
		void join(int comm1, int comm2);

		double gain(int node, int comm);
		double swap_gain(int node, int comm1, int comm2);
		double join_gain(int comm1, int comm2);
		
		void init_cp();
		void make_label_pair();
		void print_basic();
		void print_partial_sums();		
};

///turn the couble index to 
template <typename T>
void CPModularity<T>::make_label_pair(){ 
	
	unordered_map<int, int> tmp;
	node_to_cp.clear();
	for(auto &i : node_to_comm){
		tmp[i.first] = i.second/2;
		node_to_cp[i.first] = i.second%2;
	}
	node_to_comm = tmp;
	
}






template <typename T>
CPModularity<T>::CPModularity(T& _G) : Quality<T>(_G){
	quality_type = "cp_modularity";	
	node_to_cp.clear();
	for(auto &i : G.get_nbrs){ node_to_cp[ i.first ] = 1; } 
	init_cp();
	CPModularity<T>::init();
}

template <typename T>
CPModularity<T>::CPModularity(T& _G, unordered_map<int, int> &labels) : Quality<T>(_G, labels){
	quality_type = "cp_modularity";
	node_to_cp.clear();
	for(auto &i : G.get_nbrs){ node_to_cp[ i.first ] = 1; } 
	init_cp();
	CPModularity<T>::init();
}

template <typename T>
CPModularity<T>::CPModularity(T& _G, unordered_map<int, int> &labels, unordered_map<int, int> &cp_labels) : Quality<T>(_G, labels){
	quality_type = "cp_modularity";
	node_to_cp = cp_labels;
	init_cp();
	CPModularity<T>::init();
}

template <typename T>
void CPModularity<T>::set_norm(){
	norm = 1;
}

///output basic info
template <typename T>
void CPModularity<T>::print_basic(){
	Quality<T>::print_basic();
	for(auto &i : node_to_comm){
		cout << i.first << " is " << ((i.second%2)?"core":"periphery") << endl;
	}
}

///print partial sums
template <typename T>
void CPModularity<T>::print_partial_sums(){
	for(auto &i : in){ cout << " in[" << i.first << "] = " << i.second << " inc[" << i.first << "] = " << inc[i.first] << " tot[" << i.first << "] = " << tot[i.first] << " totc[" << i.first << "] = " << totc[i.first] << endl; }
}


///Set up community maps from a labelling
template <typename T>
void CPModularity<T>::init_cp(){ 
	
	unordered_map<int, int> tmp;
	for(auto &i : node_to_comm){
		tmp[i.first] = 2*i.second + node_to_cp[i.first];
	}
	node_to_comm = tmp;
	comm_to_node.clear();
	for(auto &i : node_to_comm){
		if(comm_to_node.find(i.second) == comm_to_node.end()){
			comm_to_node[i.second] = unordered_set<int>{ i.first };
		} else {
			comm_to_node[i.second].insert( i.first );
		}
	}
	num_communities = comm_to_node.size();
}

///set up in and tot vectors
template <typename T>
void CPModularity<T>::init(){
	

	in.clear();
	inc.clear();
	tot.clear();
	totc.clear();

	for(auto &i : comm_to_node){ in[i.first/2] = 0; inc[i.first/2] = 0; tot[i.first/2] = 0; totc[i.first/2] = 0; }
	
	
	//loop over edges
	for(auto &i: G.get_nbrs){
		int xi = node_to_comm[i.first]%2;
		int ci = node_to_comm[i.first]/2;
		for(auto &j: i.second){
			int xj = node_to_comm[j.first]%2;
			int cj = node_to_comm[j.first]/2;
			if( ci == cj ){
				in[ ci ] += (xj)*j.second;
				inc[ ci ] += (xi*xj)*j.second;
			} 
		}
	}

	//loop over nodes
	for(auto &i : G.degrees){ 
		tot[ node_to_comm[i.first]/2 ] += i.second; 							
		totc[ node_to_comm[i.first]/2 ] += i.second * (node_to_comm[i.first]%2); 	
	}
	
}

template <typename T>
void CPModularity<T>::eval(){
	
	double q = 0;
	for(auto &i : in){			
		q += 2*i.second - inc[i.first] - gamma*totc[i.first]*(2*tot[i.first] - totc[i.first])/G.num_links;	
		//cout << "q += 2*" << i.second << "-" << inc[i.first] << " - " << 	totc[i.first] << "*(2*" << tot[i.first] << "-" << totc[i.first] << ")" << endl;
	}
	val = min*q/G.num_links;

}

///Calculate how many links	go from node n to the different communities
template <typename T>
void CPModularity<T>::links_to_comms(int n){
	dc.clear();
	dcc.clear();
	for(auto &i : in){ 
		dc[ i.first] = 0; 
		dcc[ i.first] = 0; 
	} 
	//links where n is core. where n is not core, these are 0
	for(auto &i : G.get_nbrs[n]){
		if(i.first != n){ //skip self loops
			int ci = node_to_comm[i.first]/2;
			int xi = node_to_comm[i.first]%2;
			dc[ ci ]  +=  i.second;
			dcc[ ci ]  +=  i.second * xi;
		} 
	}
 
}

template <typename T>
void CPModularity<T>::remove(int node, int comm){
	int c = comm/2;
	int x = comm%2;


    in[c] -= x*dc[c] + x*G.selfloop(node); 
    inc[c] -= x*dcc[c] + x*G.selfloop(node); 
	tot[c] -= G.degrees[node];
	totc[c] -= x*G.degrees[node];
		
	update_nodemaps(node, comm, -1);
	node_to_cp[node] = 1; //node is in own community as core

}

template <typename T>
void CPModularity<T>::insert(int node, int comm){

	int c = comm/2;
	int x = comm%2;
	
    in[c] += x*dc[c] + x*G.selfloop(node); 
    inc[c] += x*dcc[c] + x*G.selfloop(node); 
	tot[c] += G.degrees[node];
	totc[c] += x*G.degrees[node];
			
	update_nodemaps(node, -1, comm);
	node_to_cp[node] = x; 

}

template <typename T>
double CPModularity<T>::gain(int node, int comm) {

	int c = comm/2;
	int x = comm%2;


	return min*norm*(   (x*(2*dc[c] - dcc[c]) - G.selfloop(node)*(1-x)) - 
					gamma*( (1-x)*G.degrees[node]*(2*totc[c] -  G.degrees[node]) + x*tot[c]*2*G.degrees[node] )/G.num_links   );

}

template <typename T>
double CPModularity<T>::swap_gain(int node, int comm1, int comm2) {

	int c1 = comm1/2;
	int x1 = comm1%2;
	int c2 = comm2/2;
	int x2 = comm2%2;
	
	return min*norm*(
			2*(x2*dc[c2] - x1*dc[c1]) - (x2*dcc[c2] - x1*dcc[c1]) + G.selfloop(node)*(x2 - x1)
			-gamma*(
			2*( totc[c2]* (1-x2) - totc[c1]* (1-x1) ) + 2 *(tot[c2] *x2 - tot[c1] *x1) + G.degrees[node]*( x2*(2-x2) + x1*(2-x1) ) 
			)*G.degrees[node]/G.num_links
			);
	
}

template <typename T>
void CPModularity<T>::links_between(int comm1, int comm2){
	v12 = 0;
	vc1 = 0;
	
	int c2 = comm2/2;
	int x2 = comm2%2;

	for(auto &i : comm_to_node[comm1]){
		for(auto &j : G.get_nbrs[i]){
			int cj = node_to_comm[j.first]/2;
			int xj = node_to_comm[j.first]%2;
			if(cj == c2){
				v12 += j.second * xj;
			}
			if(node_to_comm[j.first] == comm1){
				vc1 += j.second;
			}
		}
	}
 
}



template <typename T>
void CPModularity<T>::join(int comm1, int comm2){
 
	int c1 = comm1/2;
	int x1 = comm1%2;
	int c2 = comm2/2;
	int x2 = comm2%2;
	

	in[c2]  += v12 + x2*vc1;
	in[c1] = (1-x1)*(in[c1]);
	inc[c2]  += x2*(v12 + x2*vc1);
	inc[c1] = (1-x1)*(in[c1]);
	

	double d12 = (tot[c1]-totc[c1]) + x1*(2*totc[c1] - tot[c1]);
	tot[c2] += d12;
	totc[c2] += x2*d12;
	totc[c1] -= x1*totc[c1];
	tot[c1] -= d12;


	for(auto &i : comm_to_node[comm1]){ 
		node_to_comm[ i ] = comm2; 
		comm_to_node[comm2].insert( i ); 
	}
	comm_to_node.erase(comm1);	
	--num_communities;  
  
}


template <typename T>
double CPModularity<T>::join_gain(int comm1, int comm2) {
// Q = (1/2E) * sum_c ( 2*inp [c] - inc [c]  - totc[ c ]*( 2*tot[c] - totc[c] ) /2E  )

	int c1 = comm1/2;
	int x1 = comm1%2;
	int c2 = comm2/2;
	int x2 = comm2%2;
		
	double in2 = in[c2] + v12 + x2*vc1;
	double in1 = (1-x1)*(in[c1]);
	double inc2 = inc[c2] + x2*(v12 + x2*vc1);
	double inc1 = (1-x1)*(in[c1]);
	
	double in_part = (2*(in1+in2) - (inc2+inc1)) - (2*(in[c2]+in[c1]) - (inc[c2]+inc[c1]));
	
 	double d12 = (tot[c1]-totc[c1]) + x1*(2*totc[c1] - tot[c1]);
	double tot2 = tot[c2] + d12;
	double totc2 = totc[c2] + x2*d12;
	double totc1 = totc[c1] - x1*totc[c1];
	double tot1 = tot[c1] - d12;

	double tot_part = (totc2*( 2*tot2 - totc2 ) + totc1*( 2*tot1 - totc1 )) - (totc[ c2 ]*( 2*tot[c2] - totc[c2] ) + totc[ c1 ]*( 2*tot[c1] - totc[c1] ));
	
	return min*norm*( in_part - gamma*tot_part/G.num_links );
 
 
}
#endif
