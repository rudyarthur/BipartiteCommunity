#ifndef ASYNCLABELPROP_H
#define ASYNCLABELPROP_H

#include <unordered_set>
#include <set>
#include <algorithm>

#include "labelprop.h"


using namespace std; 


template <typename T>
class AsyncLabelProp : public LabelProp<T> {
	public:
		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::val;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::labels;
		//functions
		using LabelProp<T>::calc_high_freq_labels; 
		using Optimiser<T>::fixupQ; 
		using Optimiser<T>::shuffle; 
		using LabelProp<T>::use_ew; 
		using LabelProp<T>::passes; 
		
		AsyncLabelProp(T& Q, bool _use_ew=true);
		
		double optimise();
};

template <typename T>
AsyncLabelProp<T>::AsyncLabelProp(T& _Q, bool _use_ew) : LabelProp<T>(_Q, _use_ew) {
	optimiser_type = "async_labelprop";
}

///pretty dumb algorithm for label prop
template <typename T>
double AsyncLabelProp<T>::optimise() {
	
	passes = 0;
	int moves;
	do{
		
		vector<int> node_list; for(auto &i : Q.node_to_comm){ node_list.push_back(i.first); } shuffle(node_list); //random order
		moves = 0;
		for(auto &node : node_list){
			pair<int, int> high = calc_high_freq_labels(node);
			if( !high.second ){ 
				Q.node_to_comm[node] = high.first; //could be a random choice but this ought to do!
				++moves; 
			}
		}
		++passes;
		
	} while (moves > 0);
	
	//fix the community lists
	fixupQ();
	labels = Q.node_to_comm;
	Q.eval();
	
	return 0;
}


#endif
