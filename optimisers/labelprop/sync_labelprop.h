#ifndef SYNCLABELPROP_H
#define SYNCLABELPROP_H

#include <unordered_set>
#include <set>
#include <algorithm>

#include "labelprop.h"

using namespace std; 

template <typename T>
class SyncLabelProp : public LabelProp<T> {
	public:
		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::val;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::labels;
		using LabelProp<T>::high_freq_labels; 
		//functions
		using Optimiser<T>::fixupQ;		
		using LabelProp<T>::calc_high_freq_labels; 
		using LabelProp<T>::use_ew; 
		using LabelProp<T>::passes; 

		
		//unordered_map<int, unordered_set<int> > colouring;
		
		SyncLabelProp(T& Q, bool _use_ew=false);
		
		void colour();
		double optimise();
};

template <typename T>
SyncLabelProp<T>::SyncLabelProp(T& _Q, bool _use_ew) : LabelProp<T>(_Q, _use_ew) {
	optimiser_type = "sync_labelprop";
}

///Synchronous label propagation. All labels updated at same time. 
///This is the same algorithm networkx uses
template <typename T>
double SyncLabelProp<T>::optimise() {
	colour();
	
	bool stable = true;
	passes = 0;
	
	do{
		
		//for every colour
		for(auto &i : Q.comm_to_node){ //this now contains the colour groups
			//for every node
			for(auto &node : i.second){
				pair<int, bool> high = calc_high_freq_labels(node);
				//label the node by its neighbours
				if( !high.second ){ 
					Q.node_to_comm[node] = high.first; 
				} 
			}
		}
		
		//for every node
		for(int i=0; i<num_nodes; ++i){
			pair<int, bool> high = calc_high_freq_labels(i);			
			//check if there's a better label
			if( !high.second ){
				stable = false; //if so keep label propagating
				break;
			}
		}
		++passes;
		
	} while (!stable);
	//fix the community lists
	fixupQ();
	labels = Q.node_to_comm;
	Q.eval();

	return 0;
}

///Doing this to have an order to search the network where updating one node label doesn't instantly invalidate the previous ones
template <typename T>
void SyncLabelProp<T>::colour() {
	//colour all the nodes a different colour.
	//Simple but probably not the best possible implementation
	
	//define an order to go through the nodes (by degree), init all labels to -1
	Q.node_to_comm.clear(); 
	vector<int> order(num_nodes); 
	for(int i=0;i<num_nodes;++i){ order[i] = i; Q.node_to_comm[i] = -1; } 
	sort (order.begin(), order.end(), [&] (auto& a, auto& b)->bool{ return (Q.G.degrees[a] > Q.G.degrees[b]); }  );

	
	int max_col = -1;
    for(int i=0;i<num_nodes;++i){
		int node = order[i];
		
		unordered_set<int> neighbour_colours;
		for(auto &j : Q.G.get_nbrs[node]){ neighbour_colours.insert( Q.node_to_comm[ j.first ] ); } 
		
		//find the smallest 'colour' not yet seen among the neighbours
		int c = 0;
		for(; ; ++c){
			if(neighbour_colours.find(c) == neighbour_colours.end()){break;}
		}
		//is this the largest colour so fat?
		if( c>max_col ){ max_col=c; }
		//relabel node
		Q.node_to_comm[node] = c;
	}


	Q.comm_to_node.clear(); //this holds the colors
	for(auto &i : Q.node_to_comm){
		if(Q.comm_to_node.find(i.second) == Q.comm_to_node.end()){
			Q.comm_to_node[i.second] = unordered_set<int>{ i.first };
		} else {
			Q.comm_to_node[i.second].insert( i.first );
		}
	}
	for(auto &i : Q.G.get_nbrs){ Q.node_to_comm[i.first] = i.first; } //these are the initial labels



}
       
        

#endif
