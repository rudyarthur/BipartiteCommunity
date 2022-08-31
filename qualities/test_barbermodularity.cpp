#include <iostream>
#include "../networks/network_utils.h"
#include "../networks/weighted_network.h"
#include "quality.h"
#include "barber_modularity.h"
#include <unordered_map>


using namespace std; 


int main(int argc, char **argv){
	
	cout << ">> Empty network" << endl;
	BipartiteNetwork G;
	G.print_basic();
	
	cout << ">> read network from file" << endl;
	print_header("../test_data/bipartite1.txt");
	G.read("../test_data/bipartite1.txt");
	
	//checking the basic louvain methods
	{
		unordered_map<int,int> labels{ {0,1},{1,1},{2,1},{3,2},{4,3},{5,3},{6,1},{7,1},{8,2},{9,3} };
		BarberModularity Q(G,labels);	
		Q.eval();
		Q.print_basic();
		Q.print_partial_sums();

		cout << ">> Calculate number of links to each community" << endl;
		print_header("../test_data/bipartite1.txt");
		for(auto &i : Q.node_to_comm){
			Q.links_to_comms( i.first );
			for(auto &j : Q.dc){
				cout << "node " << i.first << " --[" << j.second << "]-- community " << j.first << endl;
			}
		}
	
	

		cout << ">> Remove node 0 from community 1" << endl;
		Q.links_to_comms(0);
		
		//remove
		Q.remove(0,1);	
		//quick fix to calculate checks
		Q.in[-1] = 0;
		Q.tot[-1] = Q.G.degrees[0];
		Q.totr[-1] = 0;
		Q.eval();
		double Qbefore = Q.val;
		//Q.print_basic();
		//Q.print_partial_sums();
		
			
		double dQ = Q.gain(0,1)/Q.F;	
		cout << "gain " << dQ << endl;
		cout << ">> Insert node 0 into community 1" << endl;
		
		Q.insert(0,1);
		//quick fix to calculate checks
		Q.in.erase(-1);
		Q.tot.erase(-1);
		Q.totr.erase(-1);
		Q.eval();
		double Qafter = Q.val;
		//Q.print_basic();
		//Q.print_partial_sums();
			
		cout << "dQ = " << dQ << " = " << Qafter - Qbefore << " = " << Qafter << " - " << Qbefore << endl;
	}

	//direct swaps
	{
		
		unordered_map<int,int> labels2{ {0,1},{1,1},{2,1},{3,2},{4,3},{5,3},{6,1},{7,1},{8,2},{9,3} };
		BarberModularity Q2(G,labels2);	
		//Q2.print_partial_sums();
		Q2.eval();
		double Qtarget = Q2.val;
		
		unordered_map<int,int> labels{ {0,3},{1,1},{2,1},{3,2},{4,3},{5,3},{6,1},{7,1},{8,2},{9,3} };
		BarberModularity Q1(G,labels);	
		Q1.eval();
		double Qbefore = Q1.val;

		cout << ">> Swapping node 0 from 3 into 1" << endl;
		Q1.links_to_comms(0);
		double dQ = Q1.swap_gain(0, 3, 1)/Q1.F;
		Q1.swap(0, 3, 1);
		//Q1.print_partial_sums();
		Q1.eval();
		double Qafter = Q1.val;


		cout << Qtarget << " = " << Qafter << endl;
		cout << "dQ = " << Qafter << " - " << Qbefore << " = " << Qafter - Qbefore << " = " << dQ << endl;
	}
	
	//merges
	{
		
		unordered_map<int,int> labels2{ {0,1},{1,1},{2,1},{3,3},{4,3},{5,3},{6,1},{7,1},{8,3},{9,3} };
		BarberModularity Q2(G,labels2);	
		//Q2.print_partial_sums();
		Q2.eval();
		double Qtarget = Q2.val;
		
		unordered_map<int,int> labels{ {0,1},{1,1},{2,1},{3,2},{4,3},{5,3},{6,1},{7,1},{8,2},{9,3} };
		BarberModularity Q1(G,labels);	
		Q1.eval();
		double Qbefore = Q1.val;

		cout << ">> merging community 3 and 2" << endl;
		Q1.links_between(2,3);
		double dQ = Q1.join_gain(2, 3)/Q1.F;
		
		Q1.join(2,3);
		Q1.eval();
		//Q1.print_basic();
		//Q1.print_partial_sums();
		double Qafter = Q1.val;
		
		cout << Qtarget << " = " << Qafter << endl;
		cout << "dQ = " << Qafter << " - " << Qbefore << " = " << Qafter - Qbefore << " = " << dQ << endl;


	}
	return 0;
}
