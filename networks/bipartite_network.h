#ifndef BIPARTITENETWORK_H
#define BIPARTITENETWORK_H

#include "network.h"
#include <queue>

using namespace std; 

class BipartiteNetwork : public Network {
	public:
		

		int B; //indexing will be left node:(0,1, ..., B-1) right nodes:(B, B+1, ..., num_nodes-1)
		map<int, int> colour;	//colour vector for bipartite check
					
		BipartiteNetwork():Network(){
			network_type="bipartite_network";
			B = 0;
		}		
		BipartiteNetwork(const char* filename);
		BipartiteNetwork(adjacency_map &lk);
		BipartiteNetwork(adjacency_map &lk, int node, int init_colour);
		BipartiteNetwork(const Network &G );
		void read(const char* filename);

		void print_basic(bool use_original_ids=false);
		
		//bool check_bipartite(map< int, vector<link> > &lk);
		bool is_bipartite(int init_node, int init_colour);
		void compute_indices();
		void reindex();
		BipartiteNetwork induced_graph(unordered_map<int, int> &label, int side, bool unweighted=false);

		
};

void BipartiteNetwork::read(const char* filename){
	Network::read(filename);	
	if(not is_bipartite( get_nbrs.begin()->first, 0 )){ cerr << filename << " does not contain a bipartite network!" << endl; }
}	
BipartiteNetwork::BipartiteNetwork(const char* filename){
	network_type="bipartite_network";	
	BipartiteNetwork::read(filename);
}

BipartiteNetwork::BipartiteNetwork(adjacency_map &mp) : Network(mp){
	network_type="bipartite_network";	
	if(not is_bipartite( get_nbrs.begin()->first, 0 )){ cerr << "mp does not contain a bipartite network!" << endl; }	
}
BipartiteNetwork::BipartiteNetwork(adjacency_map &mp, int node, int init_colour) : Network(mp){
	network_type="bipartite_network";	
	if(not is_bipartite( node, init_colour )){ cerr << "mp does not contain a bipartite network!" << endl; }	
}

BipartiteNetwork::BipartiteNetwork(const Network &G) : Network(G){
	network_type="bipartite_network";	
	if(not is_bipartite( get_nbrs.begin()->first, 0 )){ cerr << "G does not contain a bipartite network!" << endl; }	
}


///print some stuff
void BipartiteNetwork::print_basic(bool use_original_ids){

	Network::print_basic(use_original_ids);
	for(auto &c : colour){ cout << "colour " << c.first << " " << c.second << endl; }

	cout << "Number of left nodes = " << B << endl;
	cout << "Number of right nodes = " << num_nodes - B << endl;
}


///See if a map is bipartite
bool BipartiteNetwork::is_bipartite(int init_node, int init_colour){
	
	//colouring algorithm, will not work for disconnected networks
	queue<int> q; q.push( init_node );
	colour.clear(); colour[ init_node ] = init_colour;
	
	while(!q.empty()){
		int i = q.front(); q.pop();
		for(auto &j : get_nbrs[i]){ 
			int nn = j.first; //the neighbour node
			if(colour.find(nn) != colour.end()){ //have we already coloured it?
				if(colour[i] == colour[nn]){ return false; } //if two neighbours have same color, not bipartite
			} else {
				colour[nn] = 1-colour[i]; //colour the neighbour the opposite. 1-1 = 0, 1-0 = 1
				q.push(nn);				  //add neighbour to queue
			}
		}
	}
	B = num_nodes;
	for(auto &c : colour){ B -= c.second;  }
	
	return true;

}

///compute indices of bipartite network.
///put all left nodes at [0,1,...,B-1]
///right nodes at [B, B+1, ..., num_nodes-1]
void BipartiteNetwork::compute_indices(){
	

	new_to_old.resize(num_nodes);	
	B = 0;	
	for(auto &c : colour){
		if( c.second == 0 ){ new_to_old[B] = c.first;  old_to_new[c.first] = B; ++B; }
	}
	int C = B;
	for(auto &c : colour){
		if( c.second != 0 ){ new_to_old[C] = c.first;  old_to_new[c.first] = C; ++C; }
	}

}

///compute the node degrees & neighbour list
void BipartiteNetwork::reindex(){
	
	compute_indices(); //generate the new indices
	copy_nbrs();
	reindexed = true;
	
}


BipartiteNetwork BipartiteNetwork::induced_graph(unordered_map<int, int> &label, int side, bool unweighted){
	

	adjacency_map mp;
	int init_node=0;
	int init_colour=0;


	for(auto &i : label){	
		int a = ( colour[i.first] == side ) ? i.second : i.first;
		if(colour[i.first] == side){ init_node = i.second ; init_colour = side; }
		for(auto &j : get_nbrs[i.first]){
			int b = ( colour[j.first] == side ) ? label[j.first] : j.first;

			if(mp.find(a) == mp.end()){ 
				mp[a][b] = (unweighted) ? 1 : j.second; 
				mp[b][a] = (unweighted) ? 1 : j.second; 
			}
			else{
				if(mp[a].find(b) == mp[a].end()){ 
					mp[a][b] = (unweighted) ? 1 : j.second;  
					mp[b][a] = (unweighted) ? 1 : j.second;  
				}
				else{	
					mp[a][b] += (unweighted) ? 0 : j.second; 
					mp[b][a] += (unweighted) ? 0 : j.second;  
				}
			}
		}		
	}
	BipartiteNetwork bp = BipartiteNetwork(mp, init_node, init_colour);
	
	return bp;
}




#endif
