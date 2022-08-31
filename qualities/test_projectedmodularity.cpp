#include <iostream>
#include "../networks/network_utils.h"
#include "../networks/weighted_network.h"
#include "../networks/bipartite_network.h"
#include "../networks/projector.h"
#include "quality.h"
#include "projected_modularity.h"
#include <unordered_map>


using namespace std; 


int main(int argc, char **argv){
	
	cout << ">> Empty network" << endl;
	BipartiteNetwork G;
	G.print_basic();
	
	cout << ">> read network from file" << endl;
	print_header("../test_data/bipartite1.txt");
	G.read("../test_data/bipartite1.txt");
	
	WeightedProjector P1(G, 0, false);
	WeightedProjector P2(G, 1, false);
	WeightedProjector P3(G, 0, true);
	WeightedProjector P4(G, 1, true);
	vector<WeightedProjector> Proj = {P1, P2, P3, P4};
	//for(unsigned i=0; i<Proj.size(); ++i){ Proj[i].project(); }
	
	vector< unordered_map<int,int> > labels = {
		{ {0,1},{1,1},{2,1},{3,2},{4,2},{5,2} },
		{ {6,1},{7,1},{8,2},{9,2} },
		{ {0,1},{1,1},{2,1},{3,2},{4,2},{5,2} },
		{ {6,1},{7,1},{8,2},{9,2} }
	};
	
	vector< unordered_map<int,int> > slabels = {
		{ {0,3},{1,1},{2,1},{3,2},{4,2},{5,2} },
		{ {6,1},{7,3},{8,2},{9,2} },
		{ {0,3},{1,1},{2,1},{3,2},{4,2},{5,2} },
		{ {6,1},{7,3},{8,2},{9,2} }
	};


			
	vector< int > move = {0,7,0,7};

	//checking the basic louvain methods
	for(unsigned i=0; i<4; ++i){
		
		cout << "Projector " << i+1 << endl;
		ProjectedModularity<WeightedProjector> Q(Proj[i],labels[i]);
		//Q.G.print_basic();	
		//Q.eval();
		//Q.print_basic();
		//Q.print_partial_sums();

		cout << ">> Remove node " <<  move[i] << " from community 1" << endl;
		Q.links_to_comms(move[i]);
	
		//remove
		Q.remove(move[i],1);	
		//quick fix to calculate checks
		Q.in[-1] = Q.G.selfloop( move[i] );
		Q.tot[-1] = Proj[i].G.degrees[ move[i] ];
		Q.eval();
		double Qbefore = Q.val;
		//Q.print_basic();
		//Q.print_partial_sums();
	
		
		double dQ = Q.gain(move[i],1)*(2/Proj[i].GP.num_links) ;	
		cout << ">> Insert node " << move[i] << " into community 1" << endl;
	
		Q.insert(move[i],1);
		//quick fix to calculate checks
		Q.in.erase(-1);
		Q.tot.erase(-1);
		Q.eval();
		double Qafter = Q.val;
		//Q.print_basic();
		//Q.print_partial_sums();
		
		cout << "dQ = " << dQ << " = " << Qafter - Qbefore << " = " << Qafter << " - " << Qbefore << endl;
	}
	
	//direct swaps
	for(unsigned i=0; i<Proj.size(); ++i){
		
		cout << "Projector " << i+1 << endl;
		ProjectedModularity<WeightedProjector> Q2(Proj[i],labels[i]);
		Q2.eval();
		double Qtarget = Q2.val;
		
		ProjectedModularity<WeightedProjector> Q1(Proj[i],slabels[i]);	
		Q1.eval();
		double Qbefore = Q1.val;

		cout << ">> Swapping node " << move[i] << " from 3 into 1" << endl;
		Q1.links_to_comms(move[i]);
		double dQ = Q1.swap_gain(move[i], 3, 1)*(2/Proj[i].GP.num_links) ;
		Q1.swap(move[i], 3, 1);
		Q1.eval();
		double Qafter = Q1.val;


		cout << Qtarget << " = " << Qafter << endl;
		cout << "dQ = " << Qafter << " - " << Qbefore << " = " << Qafter - Qbefore << " = " << dQ << endl;
	}
	
	//merges
	for(unsigned i=0; i<Proj.size(); ++i){
		
		cout << "Projector " << i+1 << endl;
		ProjectedModularity<WeightedProjector> Q2(Proj[i],labels[i]);
		Q2.eval();
		double Qtarget = Q2.val;
		

		ProjectedModularity<WeightedProjector> Q1(Proj[i],slabels[i]);	
		Q1.eval();
		double Qbefore = Q1.val;
		//Q1.print_basic();
		
		cout << ">> merging community 3 and 1" << endl;
		Q1.links_between(1,3);
		//cout << "v12 = " << v12 << endl; 
		double dQ = Q1.join_gain(1, 3)*(2/Proj[i].GP.num_links) ;
		
		Q1.join(1,3);
		Q1.eval();
		//Q1.print_basic();
		//Q1.print_partial_sums();
		double Qafter = Q1.val;
		
		cout << Qtarget << " = " << Qafter << endl;
		cout << "dQ = " << Qafter << " - " << Qbefore << " = " << Qafter - Qbefore << " = " << dQ << endl;


	}
	
	return 0;
}
