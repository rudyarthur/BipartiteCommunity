#ifndef NETWORK_H
#define NETWORK_H

#ifdef DEBUG 
	#define D(x) (x)
#else 
	#define D(x) do{}while(0)
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <set>

using namespace std; 

typedef pair<int, double> link;
typedef unordered_map< int, unordered_map<int, double> > adjacency_map;

///Interface for network objects
class Network {
	public:
		int num_nodes; 		
		double num_links;

		unordered_map< int, double > degrees;			///node degrees
		adjacency_map get_nbrs;							///neighbours of a node
		
		vector<int> new_to_old; 			///mapping between new and old node ids
		map<int, int> old_to_new;			///mapping between old and new node ids
		
		string network_type;				///Identifier
		bool reindexed;						///have the indices been remapped?
		
		Network(){
			num_nodes = 0;
			num_links = 0;
			reindexed = false;
			degrees.clear();
			get_nbrs.clear();
			network_type="network";
		}		
		

		Network(adjacency_map &mp);
		void set_from_map(adjacency_map &mp);
		void read(const char* filename);

		void compute_degrees();
		void copy_nbrs();
		void compute_indices();
		void reindex();

		void print_adjacency();
		void print_indexmap();
		void print_connections(bool use_original_ids=false);	
		void print_basic(bool use_original_ids=false);
		adjacency_map to_map();
	
		double selfloop(int node);
		Network induced_graph(unordered_map<int, int> &label, bool unweighted); //need to implement this for each type
		void join(int i, int j);
		
};


/**
 * combine node j with node i
 * i -> (i,j)
 * if k != i
 * 	A_{ik} = A_{ik} + A_{jk}
 * if k == i
 *  A_{ii} = A_{ii} + A_{ij}
**/
void Network::join(int i, int j){ 
	if(i != j){
		for(auto &k : get_nbrs[j]){ //k is neighbour of j
			int a = k.first;
			if(a == i){ //connections between i and j become self loops
				if( get_nbrs[i].find(a) == get_nbrs[i].end() ){ //does i already have a self loop?
					get_nbrs[i][i] = 2*k.second;		//if not add one
				} else {
					get_nbrs[i][i] += 2*k.second;		//otherwise add the ij link
				}
				get_nbrs[a].erase(j);  					
			} else if(a == j) { //self loops of j become self loops of i
				if( get_nbrs[i].find(a) == get_nbrs[i].end() ){
					get_nbrs[i][i] = k.second;
				} else {
					get_nbrs[i][i] += k.second;
				}
			} else { //neighbours of j are now neighbours of i
				if( get_nbrs[i].find(a) == get_nbrs[i].end() ){ //was k already a neighbour? 
					get_nbrs[i][a] = k.second;	//if not add it
					get_nbrs[a][i] = k.second;
				} else {
					get_nbrs[i][a] += k.second; //otherwise increase the edge weight
					get_nbrs[a][i] += k.second;
				}
				get_nbrs[a].erase(j);  
			}

		}
		//final purge of j
		degrees[i] += degrees[j];
		get_nbrs.erase(j);  
		degrees.erase(j);
		--num_nodes;
	} else {
		cerr << "Network::join - trying to merge a node with itself" << endl;
	}
}

Network::Network(adjacency_map &mp){
	network_type="network";		
	set_from_map(mp);
}

///Set up a network object from a map
void Network::set_from_map(adjacency_map &mp){
	
	get_nbrs = mp;
	num_nodes = mp.size();
	num_links = 0; 
	for(auto it = mp.begin(); it != mp.end(); ++it){ 
		for (auto &j : it->second) {
			num_links += j.second ;
		} 
	}
	compute_degrees();
}

///return a (new) map representation of the Network
adjacency_map Network::to_map(){
	
	adjacency_map mp;

	for(auto &i : get_nbrs){
		unordered_map<int, double> tmp;
		for( auto &j : i.second ){ tmp[j.first] = j.second; }
		mp[i.first] = tmp;
	}
	
	return mp;
	
}

///compute the node degrees
void Network::compute_degrees(){

	degrees.clear();
	for(auto &i: get_nbrs){
		double d = 0;
		for( auto &j : i.second ){ d += j.second; }
		degrees[i.first] = d;
	}
	
}

///Ensure node labels go from 0 to num_nodes-1. Keep a map to the input labels
void Network::compute_indices(){

	new_to_old.resize(0);
	old_to_new.clear();
	
	set<int> keys;
	for (auto& it : get_nbrs) keys.insert(it.first); //node names
	for( auto& it : keys){
		old_to_new[it] = new_to_old.size(); //new indices are from [0, keys.size()-1]
		new_to_old.push_back(it);
	}

}

void Network::copy_nbrs(){

	//recompute everything
	degrees.clear();	
	adjacency_map tmp;
	
	for(unsigned i=0; i<get_nbrs.size(); ++i){ 

		double d = 0;
		for( auto &j : get_nbrs[ new_to_old[i] ] ){ d += j.second; }
		degrees[i] = d;
		
		for( auto &j : get_nbrs[ new_to_old[i] ] ){
			tmp[i][ old_to_new[ j.first ] ] = j.second;
		} 
		
	}
	
	get_nbrs = tmp;

}

///compute the node degrees & neighbour list
void Network::reindex(){
	
	compute_indices(); //generate the new indices
	copy_nbrs();
	reindexed = true;
	
}
	

/**
 * Read a network in from a file
 * Expected format
 * weighted = false
 * 0 1
 * 0 2
 * 1 2 ...
 * weighted = true
 * 0 1 1
 * 0 2 1
 * 1 2 3 ...
 * Edges are assumed to be bidirectional!
 * **/
void Network::read(const char* filename){
	
	ifstream infile(filename);
	if (infile.is_open() != true) {
		cerr << "The file " << filename << " does not exist" << endl;
		exit(1);
	}
  	
	string line;

	num_links = 0;		
	while(getline(infile, line)) {
		//skip comments and blank lines: %, #, // allowed comment styles
		if(line.empty() || (line.find("%") == 0) || (line.find("#") == 0) || (line.find("//") == 0) ) {continue;} 
		
		//read the other lines
		istringstream iss(line);
		string token;
		int num_tokens = 0;

		int a, b;
		vector<int> n(2);
		double w=1;
		
		while(getline(iss, token, ' ')){
			if(num_tokens < 2){
				n[num_tokens] = atoi(token.c_str());
			} else {
				w = atof(token.c_str()); //last thing written on the line is the weight, otherwise default to 1
			}
			++num_tokens;
		}
		
		a = n[0]; b = n[1];
		get_nbrs[a][b] = w;
		if(a != b){ get_nbrs[b][a] = w; }
		else{ get_nbrs[a][a] += w; }	//this is basically a convention...

		num_links+=2*w;					//...which makes this work
			
	}
	infile.close();
	num_nodes = get_nbrs.size();
  	compute_degrees();

  	
}




///return the self edge weight at n. 
double Network::selfloop(int n){

	auto it = get_nbrs[n].find(n);
	if(it != get_nbrs[n].end()){
		return it->second;
	}
	return 0;
}

///Use a node labelling to compute an induced network	
Network Network::induced_graph(unordered_map<int, int> &label, bool unweighted){
	cerr << "induced_graph should not be called from Base class Network" << endl;
	exit(1);
	return Network();
}

///print adjacency matrix with zeros
void Network::print_adjacency(){

	set<int> keys;
	for (auto& it : get_nbrs){ keys.insert(it.first); }
	cout << " |"; for( auto& it : keys){ cout << it << " "; } cout << endl;
	cout << " |"; for( auto& it : keys){ cout << "--"; } cout << endl;
	for( auto& i : keys){ 
		cout << i << "|"; 
		for( auto& j : keys){ 
			cout << ( ( get_nbrs[i].find(j) == get_nbrs[i].end() ) ? 0 : get_nbrs[i][j] ) << " ";
		} cout << endl;
	}
}

///print the index map
void Network::print_indexmap(){
	
	cout << "Network ids have" << ((reindexed) ? "" : "NOT") << " been re-mapped" << endl;
	
	for(unsigned i = 0; i < new_to_old.size(); ++i){ 
		cout << "old_idx " << new_to_old[i] << " --> " << i << endl; 
	}

}

///print connections
void Network::print_connections(bool use_original_ids){
	for(auto &ii : get_nbrs){ 
		int a = (use_original_ids) ? new_to_old[ii.first] : ii.first;

		for( auto &j : ii.second ){
			int b = (use_original_ids) ? new_to_old[j.first] : j.first;
			cout << a << "---[" << j.second << "]---" << b << endl;
		}
	}
}

void Network::print_basic(bool use_original_ids){
	cout << "This is a " << network_type << endl;
	cout << num_nodes << " nodes" << endl;
	cout << num_links << " total edge weight" << endl;

	print_connections(use_original_ids);

	for(auto &ii : get_nbrs){ 
		int a = (use_original_ids) ? new_to_old[ii.first] : ii.first;
		cout << "neighbours of " << a << " ( ";
		for( auto &j : ii.second ){
			int b = (use_original_ids) ? new_to_old[j.first] : j.first;
			cout << b << " ";
		} cout << ")" << endl;
	}
	double cum_degree = 0;

	for(auto &i: degrees){ 
		int a = (use_original_ids) ? new_to_old[i.first] : i.first;
		cum_degree += i.second; 
		cout << "degree " << a << ": " << i.second << " " << cum_degree << endl; 
	}
}


#endif
