#ifndef DUALPROJECTION_H
#define DUALPROJECTION_H


#include "../optimiser.h"
#include "aggregate.h"
#include "projected_louvain.h"

using namespace std; 

template <typename T>
class DualProjection : public Optimiser< BarberModularity > {
	public:
		//data
		//using Optimiser<ProjectedModularity<T> >::Q; 
		//using Optimiser<ProjectedModularity<T> >::optimiser_type;
		//using Optimiser<ProjectedModularity<T> >::fixupQ;
		bool max_agg;
		bool include_loops;

		DualProjection( BarberModularity& _Q);		
		double optimise();
};

template <typename T>
DualProjection<T>::DualProjection( BarberModularity& _Q) : Optimiser< BarberModularity >(_Q) {
		optimiser_type = "dual_projection";
		max_agg = true;
		include_loops = false;
}		

template <typename T>
double DualProjection<T>::optimise() {
	
	
	double dq = 0;
	int level = 0;
	vector< unordered_map< int, int > > community_maps;

	for(int side = 0; side < 2; ++side){

		BipartiteNetwork B(Q.G);
		T P(B, side, include_loops);
		ProjectedModularity<T> QS(P);			
		ProjectedLouvain< ProjectedModularity<T> > L(QS);
		L.optimise();
		//cout << "on side " << side << endl;
		//L.print_basic();
		
		for(auto &i : L.labels){ Q.node_to_comm[ i.first ] = i.second; }	
		
	} 
	
	fixupQ();
	//Q.print_basic();


		
				
	Aggregate< BarberModularity > A( Q, max_agg );
	A.optimise();	
	//A.print_basic();

	labels.clear();
	for(auto &i: A.labels){ labels[ i.first ] = i.second; }
	
	return 0;
}


#endif
