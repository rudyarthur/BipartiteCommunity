#include <iostream>
#include "../networks/network_utils.h"
#include "../networks/weighted_network.h"
#include "quality.h"
#include "cp_modularity.h"
#include <unordered_map>


using namespace std; 


int main(int argc, char **argv){
	
	cout << ">> Empty network" << endl;
	WeightedNetwork G;
	
	cout << ">> read network from file" << endl;
	print_header("../test_data/graph5.txt");
	G.read("../test_data/graph5.txt");
	
	//unordered_map<int,int> labels   { {0,1},{1,1},{2,1},{3,1},{4,1},{5,2},{6,2},{7,2} };
	//unordered_map<int,int> cp_labels{ {0,1},{1,1},{2,1},{3,0},{4,0},{5,1},{6,1},{7,1} };

	unordered_map<int,int> labels   { {0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7} };
	unordered_map<int,int> cp_labels{ {0,1},{1,1},{2,1},{3,0},{4,0},{5,1},{6,1},{7,1} };
	
	int move = 1;
	int tocomm = 1;
	//checking the basic louvain methods
	{
		
		CPModularity<WeightedNetwork> Q(G,labels,cp_labels);			
		Q.eval();
		Q.print_basic();
		Q.print_partial_sums();
	
		//remove
		int comm = Q.node_to_comm[move];
		int c = Q.node_to_comm[move]/2;
		int x = Q.node_to_comm[move]%2;
		cout << ">> Remove node " <<  move << " from community " << c << " with cp " << x << endl;

		Q.links_to_comms(move);

		Q.remove(move,comm);	
		//quick fix to calculate checks
		Q.in[-1] = Q.G.selfloop( move );
		Q.inc[-1] = Q.G.selfloop( move );
		Q.tot[-1] = Q.G.degrees[ move ];
		Q.totc[-1] = Q.G.degrees[ move ];
		Q.eval();
		double Qbefore = Q.val;
		Q.print_basic();
		Q.print_partial_sums();
		
		double dQ = Q.gain(move,tocomm)*(1/Q.G.num_links) ;	
		cout << ">> Insert node " << move << " into community " << tocomm << endl;
	
		Q.insert(move,tocomm);
		//quick fix to calculate checks
		Q.in.erase(-1);
		Q.inc.erase(-1);
		Q.tot.erase(-1);
		Q.totc.erase(-1);
		Q.eval();
		double Qafter = Q.val;
		Q.print_basic();
		Q.print_partial_sums();
		
		cout << "dQ = " << dQ << " = " << Qafter - Qbefore << " = " << Qafter << " - " << Qbefore << endl;
	}
	/*
	//direct swaps
	{
		CPModularity<WeightedNetwork> Qt(G,labels,cp_labels);
		Qt.eval();
		//Qt.print_basic();
		//cout << "target" << endl;
		//Qt.print_partial_sums();
		double Qtarget = Qt.val;
		
		CPModularity<WeightedNetwork> Q1(G,slabels,cp_slabels);	
		Q1.eval();
		double Qbefore = Q1.val;
		//Q1.print_basic();
		//cout << "start" << endl;
		//Q1.print_partial_sums();
		
		int comm2 = Qt.node_to_comm[0];
		int c2 = Qt.node_to_comm[0]/2;
		int x2 = Qt.node_to_comm[0]%2;
		
		int comm1 = Q1.node_to_comm[0];
		int c1 = Q1.node_to_comm[0]/2;
		int x1 = Q1.node_to_comm[0]%2;

		cout << ">> Swapping node " <<  0 << " from " << comm1 << " to " << comm2 << endl;


		Q1.links_to_comms(0);
		double dQ = Q1.swap_gain(0, comm1, comm2)*(1/Q1.G.num_links) ;
		Q1.swap(0, comm1, comm2);
		Q1.eval();
		double Qafter = Q1.val;
		//cout << "after" << endl;
		//Q1.print_basic();
		//Q1.print_partial_sums();
		
		cout << Qtarget << " = " << Qafter << endl;
		cout << "dQ = " << Qafter << " - " << Qbefore << " = " << Qafter - Qbefore << " = " << dQ << endl;
	}
	
	//merges
	{
		
		unordered_map<int,int> labels2   { {0,1},{1,2},{2,1},{3,2},{4,2},{5,2},{6,3},{7,3} };
		unordered_map<int,int> cp_labels2{ {0,0},{1,0},{2,0},{3,1},{4,0},{5,0},{6,1},{7,0} };
	
		CPModularity<WeightedNetwork> Q2(G,labels2, cp_labels2);	
		Q2.eval();
		//cout << ">>Target" << endl;
		//Q2.print_basic();
		//Q2.print_partial_sums();
		double Qtarget = Q2.val;

		//unordered_map<int,int> labels   { {0,1},{1,1},{2,1},{3,2},{4,2},{5,2},{6,3},{7,3} };
		unordered_map<int,int> cp_labels{ {0,0},{1,1},{2,0},{3,1},{4,0},{5,0},{6,1},{7,0} };
		CPModularity<WeightedNetwork> Q1(G,labels,cp_labels);	
		
		Q1.eval();
		//cout << ">>Before" << endl;		
		//Q1.print_basic();
		//Q1.print_partial_sums();
		double Qbefore = Q1.val;

		cout << ">> merging core of 1 with periphery of 2" << endl;
		Q1.links_between(3,4);

		double dQ = Q1.join_gain(3,4)/G.num_links;		
		Q1.join(3,4);
		Q1.eval();
		//cout << ">>After" << endl;				
		//Q1.print_basic();
		//Q1.print_partial_sums();
		double Qafter = Q1.val;
		
		cout << Qtarget << " = " << Qafter << endl;
		cout << "dQ = " << Qafter << " - " << Qbefore << " = " << Qafter - Qbefore << " = " << dQ << endl;


	}*/
	
	return 0;
}
