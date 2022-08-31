#ifndef PROJECTION_H
#define PROJECTION_H

#include "network.h"
#include "weighted_network.h"
#include "bipartite_network.h"
#include <algorithm>

using namespace std; 

/**
 * q_i = sum_k B_ik
 * d_k = sum_i B_ik
 * 
 * Have to deal with some complications
 * Simplest projection is A_ij = sum_k B_ik B_jk
 * but usually we force A_ii = 0
 * so
 * A_ij (no loops) = sum_k B_ik B_jk - delta_ij sum_k B_ik B_jk = (1 - delta_ij)sum_k B_ik B_jk
 * generally
 * A_ij = (1 - s*delta_ij)sum_k B_ik B_jk
 * s = 0 -> loops on
 * s = 1 -> loops off
 * Add weight
 * A_ij = (1 - s*delta_ij)sum_k B_ik B_jk / w_k
 * 
 * 
 * 
 * 
 * 
 * weighted proj: w_k = 1
 * sum_ij A_ij = sum_ij (1 - s*delta_ij)sum_k B_ik B_jk / w_k
 * = sum_ijk B_ik B_jk - s*sum_ij delta_ij (sum_k B_ik B_jk)
 * = sum_k d_k^2 - s*sum_ik B_ik^2
 * = 2E (this is the definition of E)
 * 
 * 
 * 
 * 
 * Need to understand what happens in the configuration model when projecting
 * Prob(i -> k) = q_i d_k / F
 * Prob(k -> j) tricky if we allow j=i, then we have loops in the graph, if not, we don't
 * loops
 * Prob(k -> j) = d_k * q_j / F 
 * no loops
 * Prob(k -> j) = (d_k - 1) * q_j / (F - 1)  because one edge choice is removed.
 * generally
 * Prob(k -> j) = (d_k - s) * q_j / (F - s) 
 * 
 * 
 * 
 * P_ij = sum_k Prob(i -> k) Prob(k -> j) =  (sum_k  d_k (dk - s)/w_k )  *  (q_i q_j)/(F*(F-s))
 * 
 * 
 * D2 = (sum_k  d_k (dk - s)/w_k )
 * 
 * loops, weighted proj, s=0
 * D2 = sum_k  d_k d_k = 2E
 * 
 * no loops, weighted proj, s=1
 * D2 = sum_k d_k^2 - sum_k dk
 * D2 = sum_k d_k^2 - sum_ik B_ik
 * 2E = sum_k d_k^2 - sum_ik B_ik^2
 * equal if B_ik = 1
 * 
 * 
 * 
 * no loops, hyperbolic, s=1
 * D2 = sum_k d_k^2/(d_k-1) - sum_k dk/(d_k-1)
 * 2E = sum_k d_k^2/(d_k-1) - sum_ik B_ik^2/(d_k-1)
 * loops, hyperbolic, s=0
 * D2 = sum_k  d_k d_k/d_k = sum_k d_k = F = 2E
**/
class Projector {
	public:
		
		bool include_loops;	///should we include self loops (nightmare)
		int side;			///project on left?
				
		int inc;			///0 or 1 if we include or exclude loops
		BipartiteNetwork &G;	///The network we are projecting
		WeightedNetwork GP;		///The projection
		double D2;
		
		Projector(BipartiteNetwork &_G) : G(_G) { side = 0; include_loops = false; inc = 1-(int)include_loops;  }
		Projector(BipartiteNetwork &_G, int _side, bool _include_loops) : G(_G) { side = _side; include_loops = _include_loops; inc = 1-(int)include_loops;  }
		
		
		///Compute B.B^T
		adjacency_map compute_common(unordered_map<int, double> &weight){
			
			adjacency_map mp;
			for(auto &i : G.get_nbrs){	
				int a = i.first;
				if (G.colour[a] != side ){ continue; }	//skip indices on the other side
				for(auto &j : G.get_nbrs){	
					int b = j.first;	
					if (G.colour[b] != side ){ continue; }	//skip indices on the other side
					if (!include_loops && a == b){ continue; }	//don't look at self connections
					
					double wval = 0;
					bool w = false;
					//sum_k B_ik B_jk
					for(auto &k1 : G.get_nbrs[a]){	//neighbours are always on the other side!
					for(auto &k2 : G.get_nbrs[b]){
						if(k1.first == k2.first){ 
							wval += k1.second * k2.second / weight[k1.first]; 
							w = true;
						}
					}} 
					if( w ){
						mp[a][b] = wval; 
					}
		
				}
			}
			
			return mp;
			
		}
		
		void dsum(unordered_map<int, double> &weight){
			D2 = 0;
			for(auto &i : G.degrees){	
				int a = i.first;
				if (G.colour[a] == side ){ continue; }
				D2 += G.degrees[a]*(G.degrees[a]-inc)/weight[a]; 
			}
		}	

		virtual unordered_map<int, double> get_weights() = 0;
		
		
		WeightedNetwork project(){
			unordered_map<int, double> weight = get_weights();
			
			adjacency_map mp = compute_common(weight);
			dsum(weight);			
			GP = WeightedNetwork(mp);
			return GP;
		}

};

///Weighted projection
class WeightedProjector : public Projector {
	public:

		WeightedProjector(BipartiteNetwork &_G) : Projector(_G) { project(); }
		WeightedProjector(BipartiteNetwork &_G, int _side, bool _inc) : Projector(_G, _side, _inc) { project(); }
	

		unordered_map<int, double> get_weights(){
			unordered_map<int, double> weight;
			for(auto &i : G.degrees){	
				int a = i.first;
				if (G.colour[a] == side ){ continue; }
				weight[a] = 1; 
			}
			return weight;			
		}
};


///Newman projection
class HyperbolicProjector : public Projector {
	public:

		HyperbolicProjector(BipartiteNetwork &_G) : Projector(_G) { project(); }
		HyperbolicProjector(BipartiteNetwork &_G, int _side, bool _inc) : Projector(_G, _side, _inc) { project(); }
			

		unordered_map<int, double> get_weights(){
			unordered_map<int, double> weight;
			for(auto &i : G.degrees){	
				int a = i.first;
				if (G.colour[a] == side ){ continue; }
				weight[a] = G.degrees[a]-inc; 
			}
			return weight;
		}
};

#endif
