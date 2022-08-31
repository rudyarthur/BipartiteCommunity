#include <iostream>
#include "../networks/network_utils.h"
#include "../optimiser.h"
#include "projected_louvain.h"
#include <random>





using namespace std; 

int main(int argc, char **argv){

	vector< string > filenames = { "../../test_data/bipartite1.txt"};

	for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			BipartiteNetwork B(filename.c_str());
			WeightedProjector P(B, 1, false);
			
			ProjectedModularity<WeightedProjector> Q(P);			
			//Q.G.print_basic();
			//Q.P.G.print_basic();
			
			ProjectedLouvain< ProjectedModularity<WeightedProjector> > L( Q );
			L.optimise();

			print_header(filename.c_str());
			L.print_basic();

	}

	return 0;

}
