#ifndef LABELPROP_H
#define LABELPROP_H

#include "../optimiser.h"
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <algorithm>

using namespace std; 

template<typename T>      
class LabelProp : public Optimiser<T> {
	public:
		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::val;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::passes;
		using Optimiser<T>::labels;
		
		bool use_ew;			///use edge weights
		unordered_map< int, unordered_set<int> > high_freq_labels; ///for each node what label(s) appear most often among its neighbours?
		
		LabelProp(T& Q, bool _use_ew=false);
		pair<int, bool> calc_high_freq_labels(int node);		///returns a,b where a is community label with the most links to node, and b is true is node has most links to its current community
};

template <typename T>
LabelProp<T>::LabelProp(T& _Q, bool _use_ew) : Optimiser<T>(_Q) {
	optimiser_type = "labelprop";
	use_ew = _use_ew;
	passes = 0;
	high_freq_labels.clear();
}      



///get the most frequent labels of node's neighbours
template <typename T>      
pair<int, bool> LabelProp<T>::calc_high_freq_labels(int node) {
	
	
	unordered_map< int, double > freqs;
	int best_val = -1;
	for(auto &j : Q.G.get_nbrs[node]){ 
		int nbr = j.first;
		int nbr_label = Q.node_to_comm[nbr]; //node is connected to this community
		double w = (use_ew) ? j.second : 1;
		if( freqs.find(nbr_label) == freqs.end()){ freqs[ nbr_label ] = w; } //accumulate number of connections
		else { freqs[ nbr_label ] += w; }
		
		if(freqs[ nbr_label ] > best_val){ best_val = freqs[ nbr_label ]; } //community with the most connections to node
	}
	
	//get all the labels that occur as often as best_val - there can be ties
	unordered_set<int> best;
	int high_label = -1;		//the highest label that occurs most often
	int current_comm = false;	//is the current label one of the ones that occurs most often?

	for(auto &it : freqs){
		if( it.second == best_val ){ 
			best.insert(it.first); 
			if( it.first > high_label ){ high_label = it.first; }
			if( it.first == Q.node_to_comm[node] ){ current_comm = true; }
		}
	}
	
	high_freq_labels[node] = best;
	return make_pair(high_label, current_comm);
}       
        


#endif
