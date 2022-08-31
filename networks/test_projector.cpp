#include <iostream>
#include "network_utils.h"
#include "weighted_network.h"
#include "bipartite_network.h"
#include "projector.h"

using namespace std; 

int main(int argc, char **argv){

	
	BipartiteNetwork G;
	print_header("../test_data/bipartite1.txt");
	G.read("../test_data/bipartite1.txt");
	G.print_basic();
	
	WeightedProjector P_loops(G, 0, true);
	WeightedProjector P_noloops(G, 0, false);
	
	{
		WeightedNetwork GP = P_noloops.project();
		cout << "Projected without loops" << endl;
		GP.print_basic();
		cout << "D2 = " << P_noloops.D2 << endl;
	}
	{
		WeightedNetwork GP = P_loops.project();
		cout << "Projected with loops" << endl;
		GP.print_basic();
		cout << "D2 = " << P_loops.D2 << endl;
	}
		

	HyperbolicProjector HP_loops(G, 0, true);
	HyperbolicProjector HP_noloops(G, 0, false);
	
	{
		WeightedNetwork GP = HP_noloops.project();
		cout << "Newman Projected without loops" << endl;
		GP.print_basic();
		cout << "D2 = " << HP_noloops.D2 << endl;
	}
	{
		WeightedNetwork GP = HP_loops.project();
		cout << "Newman Projected with loops" << endl;
		GP.print_basic();
		cout << "D2 = " << HP_loops.D2 << endl;
	}
			
	return 0;
}
