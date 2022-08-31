#include <iostream>
#include "../networks/network_utils.h"
#include "../optimiser.h"
#include "louvain.h"
#include <random>





using namespace std; 

int main(int argc, char **argv){

	vector< string > filenames = { "../../test_data/graph1.txt" , "../../test_data/graph2.txt", "../../test_data/graph3.txt"};

	for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			WeightedNetwork W(filename.c_str());
			W.print_basic();
			Modularity<WeightedNetwork> Q(W);			
			Q.print_basic();
			
			Louvain< Modularity<WeightedNetwork> > L( Q );
			L.optimise();

			print_header(filename.c_str());
			L.print_basic();

	}

	return 0;

}
