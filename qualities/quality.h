#ifndef QUALITY_H
#define QUALITY_H

#include "../networks/network.h"
#include "../networks/weighted_network.h"
#include "../networks/bipartite_network.h"
#include "../networks/projector.h"
#include <unordered_set>
#include <unordered_map>

using namespace std; 

template <typename T>
class Quality {
	public:
		int num_communities;				//number of communities
		double val;							//value of quality fn
		double gamma;						//resolution parameter
		double norm;						//return correctly normalised gain
		double min;							//if =-1 minimise modularity
		
		unordered_map< int, int > node_to_comm;			//node to community map
		unordered_map< int, unordered_set<int> > comm_to_node;	//community to node map
				
		//the partial sums used by Louvain
		unordered_map<int, double> in;		
		unordered_map<int, double> tot;
		unordered_map<int, double> dc; //links to neighbour communities
		double v12;
		double vc1;
		
		T& G; //this must be initialised at construction
		
		string quality_type;
		
		void setup_comms(unordered_map<int, int> &labels);		
		void setup_comm_to_node();
		void setup_network(T& _G);		
		Quality(T& _G);
		Quality(T& _G, unordered_map<int, int> &labels);


		virtual void init() = 0;				
		virtual void eval() = 0;
		void set_min();
		virtual void set_norm() = 0;
		virtual void remove(int node, int comm) = 0;  //remove a node from a community
		virtual void insert(int node, int comm) = 0;	//insert a node into a communtiy
		void swap(int node, int comm1, int comm2);					//swap between communities
		virtual double gain(int node, int comm) = 0;	//(unnormalised) gain from joining an isolated node n to a communtiy
		virtual double swap_gain(int node, int comm1, int comm2) = 0;	//(unnormalised) gain from moving node from comm1 to comm2
		virtual double join_gain(int comm1, int comm2) = 0;		//(unnormalised) gain from joining comm1 and comm2
		
		
		virtual void join(int comm1, int comm2) = 0;			//join comm1 and comm2
		void qjoin(int comm1, int comm2);									//shared joining function that updates the community maps
		void update_nodemaps(int node, int comm1, int comm2);						//keep the community maps correct during remove and join
		void links_to_comms(int node);	//calculate the connections from node to all the comms		
		void links_between(int comm1, int comm2); //calculate the number of links between c1 and c2
				
		void print_basic();
		void print_partial_sums();
};

template <typename T>
void Quality<T>::set_min(){
	min = -1;
}


///Join comm1 to comm2, using the label of comm2
template <typename T>
void Quality<T>::qjoin(int comm1, int comm2){

	for(auto &i : comm_to_node[comm1]){ 
		node_to_comm[ i ] = comm2; 
		comm_to_node[comm2].insert( i ); 
	}
	comm_to_node.erase(comm1);
	in.erase(comm1);
	tot.erase(comm1);
	--num_communities;
}

/// comm1 = (.... node ... ) comm2 = (... ....) => comm1 = (... ...) comm2 = (... node ...)
template <typename T>
void Quality<T>::update_nodemaps(int node, int comm1, int comm2){
	node_to_comm[node] = comm2;
	comm_to_node[comm1].erase(node);
	comm_to_node[comm2].insert(node);
}


///Calculate how many links	go from node n to the different communities
template <typename T>
void Quality<T>::links_to_comms(int n){
	dc.clear();
	for(auto &i : comm_to_node){ dc[i.first] = 0; }
	for(auto &i : G.get_nbrs[n]){
		if(i.first != n){ //skip self loops
			dc[ node_to_comm[i.first]  ]  +=  i.second;
		} /*else if(use_loops) {
			dc[ node_to_comm[i.first]  ]  += G.selfloop(n);
		} //could maybe do this here but easier to do self loops separately*/
	}
}

///Calculate how many links	go between comm1 and comm2
template <typename T>
void Quality<T>::links_between(int comm1, int comm2){
	v12 = 0;
	for(auto &a : comm_to_node[comm1]){ 
		for(auto &j : G.get_nbrs[a]){
			int b = j.first;
			if(node_to_comm[b] == comm2){
				v12 += j.second;
			}
		}
	}
}

///output basic info
template <typename T>
void Quality<T>::print_basic(){

	cout << quality_type << " is the quality function" << endl;
	cout << num_communities << " communities" << endl;
	for(auto &i : node_to_comm){ cout << i.first << " in " << i.second << endl; }
		
	for(auto &i : comm_to_node){
		cout << "community " << i.first << " = ( ";
		for(auto &j : i.second){
			cout << j << " ";
		} cout << ")" << endl;
	}
	cout << "Current Quality = " << val << endl;
	
}

///print partial sums
template <typename T>
void Quality<T>::print_partial_sums(){
	for(auto &i : in){ cout << " in[" << i.first << "] = " << i.second << " tot[" << i.first << "] = " << tot[i.first] << endl; }
}


template <typename T>
Quality<T>::Quality(T& _G, unordered_map<int, int> &labels) : G(_G) { 
	quality_type = "quality";
	gamma = 1; 
	norm = 1;
	min = 1;	
	setup_comms(labels); 
}

///Use a network to initialise Q
template <typename T>
void Quality<T>::setup_network(T& _G){ 
	G = _G;
	unordered_map<int, int> labels; for(auto &i : G.get_nbrs){ labels[ i.first ] = i.first; }
	setup_comms( labels );	
	init();		//are we calling init twice
}

template <typename T>
Quality<T>::Quality(T& _G) : G(_G){ 

	quality_type = "quality";
	gamma = 1;
	norm = 1;
	min = 1;
	unordered_map<int, int> labels; for(auto &i : G.get_nbrs){ labels[ i.first ] = i.first; }
	setup_comms( labels );
}

///Set up community maps from a labelling
template <typename T>
void Quality<T>::setup_comms(unordered_map<int, int> &labels){ 
	
	if(labels.size() != unsigned(G.num_nodes) ){
		cerr << "Trying to use " << labels.size() << " labels with " << G.num_nodes << " node network" << endl;
		exit(1);
	}
	node_to_comm = labels;
	setup_comm_to_node();
}

template <typename T>
void Quality<T>::setup_comm_to_node(){ 
	
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

template <typename T>
void Quality<T>::swap(int node, int comm1, int comm2){

	remove( node,  comm1);
	insert( node,  comm2);
  
}

#endif
