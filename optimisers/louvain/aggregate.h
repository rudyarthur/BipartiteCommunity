#ifndef AGGREGATE_H
#define AGGREGATE_H

#include "../optimiser.h"

using namespace std; 

template <typename T>
class Aggregate : public Optimiser<T> {
	public:

		//data
		using Optimiser<T>::Q; 
		using Optimiser<T>::optimiser_type;
		using Optimiser<T>::num_nodes;
		using Optimiser<T>::num_communities;
		using Optimiser<T>::eps;
		using Optimiser<T>::passes;
		using Optimiser<T>::labels;
		
		bool max_agg; //join communities when mod change is zero
		
		unordered_map< int, unordered_map<int, double> > R; //R[c][d] = connection strength from community c to community d
		unordered_map< int, unordered_map<int, double> > Gain; //Gain[c][d] = quality gain joining community c to community d
	
		//functions
		using Optimiser<T>::fixupQ;
	
		Aggregate(T& Q, bool _ma=true);
		
		double setup_gain(int &best_comm1, int &best_comm2);
		inline void update_r(int &best_comm1, int &best_comm2);
		inline void update_gain(int &best_comm2);
		double optimise();
};


template <typename T> 
Aggregate<T>::Aggregate(T& _Q, bool _ma) : Optimiser<T>(_Q) {
	if(Q.quality_type == "cp_modularity"){
		cerr << "Aggregation is not implemented for CP Modularity" << endl;
		exit(1);
	}
	optimiser_type = "aggregate";
	max_agg = _ma;
	
	for(auto &i : Q.comm_to_node){
		for(auto &j : Q.comm_to_node){
			if(i.first <= j.first){ continue; }
			Q.links_between(i.first, j.first);
			R[i.first][j.first] = Q.v12;
			R[j.first][i.first] = Q.v12;
		} 
	}

}


template <typename T> 
double Aggregate<T>::setup_gain(int &best_comm1, int &best_comm2){

  double max_gain = 0;

  for(auto &i : Q.comm_to_node){
	int comm1 = i.first;
	for(auto &j : R[i.first]){
		
		int comm2 = j.first;
		Q.v12 = j.second; //if cp have to compute links between here for the vc1 term
		double jg = Q.join_gain(comm1, comm2);
		Gain[comm1][comm2] = jg;

		if( jg > max_gain ){
			max_gain = jg; 
			best_comm1 = comm1;
			best_comm2 = comm2;
		}

	}  
  } 
  
  return max_gain;
	
}

template <typename T> 
inline void Aggregate<T>::update_r(int &comm1, int &comm2){
	
	R[comm2].erase(comm1);
	Gain[comm2].erase(comm1);
	
	R[comm1].erase(comm2);
	Gain[comm1].erase(comm2);
		
	for(auto &it : R[comm1]){ //other eighbours of comm1 are now neighbours of comm2
		int c = it.first;	

		if( R[comm2].find(c) == R[comm2].end() ){ //comm2 was not previously linked to c
			R[comm2][c] = it.second;
			R[c][comm2] = it.second;
		} else {								  //comm2 was previously linked with c	
			R[comm2][c] += it.second;
			R[c][comm2] += it.second;				
		}	

		R[c].erase(comm1);
		Gain[c].erase(comm1); 

	}
	
	//comm1 no longer exists so remove it
	R.erase(comm1); 
	Gain.erase(comm1); 
	  
}

template <typename T> 
inline void Aggregate<T>::update_gain(int &best_comm2){ 

	for(auto &it : R[best_comm2]){
		int c = it.first;
		if(best_comm2 == c){ continue; }
		Q.v12 = it.second;   //if cp have to compute links between here for the vc1 term
		double jg = Q.join_gain(best_comm2, c);
		Gain[best_comm2][c] = jg;
		Gain[c][best_comm2] = jg;
	}  

}

template <typename T> 
double Aggregate<T>::optimise() {

	passes = 0;
	int best_comm1;
	int best_comm2;
	double max_gain = setup_gain(best_comm1, best_comm2);
	double old_val = Q.val;

	while( (max_agg) ? (max_gain >= 0) : (max_gain > eps)  ){
		
		//update community labels
		Q.join(best_comm1, best_comm2);//join comm1 and comm2, we get rid of comm1 and keep comm2!
		++passes;

		update_r(best_comm1, best_comm2);
		update_gain(best_comm2);
		
		max_gain = 0;
		best_comm1 = -1;
		best_comm2 = -1;
		for(auto &i : R){
			int comm1 = i.first;
			for(auto &j : i.second){
				int comm2 = j.first;
				if( comm1 < comm2 ){ //symmetry
					double jg = Gain[comm1][comm2];
					if( (max_agg) ? (jg >= max_gain) : (jg > max_gain) ){
						max_gain = jg; 
						best_comm1 = comm1;
						best_comm2 = comm2;
					}
				}
			}		  
		}
		if( best_comm1 == best_comm2 ){ break; }
		
	}
	fixupQ();
  	labels = Q.node_to_comm;
	
    return Q.val-old_val;
}


#endif
