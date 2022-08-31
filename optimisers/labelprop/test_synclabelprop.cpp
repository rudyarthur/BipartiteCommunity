#include <iostream>
#include "../../networks/network_utils.h"
#include "../../networks/weighted_network.h"
#include "../../qualities/quality.h"
#include "../../qualities/modularity.h"
#include "../optimiser.h"
#include "labelprop.h"
#include "sync_labelprop.h"


using namespace std; 

int main(int argc, char **argv){
	

	vector< string > filenames = { "../../test_data/graph1.txt", "../../test_data/graph2.txt", "../../test_data/graph3.txt"};

	for(auto &filename : filenames){
		cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
		cout << "\\CHECKING " << filename << endl;
		cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
		WeightedNetwork W(filename.c_str());
		Modularity<WeightedNetwork> Q(W);			
		
		SyncLabelProp< Modularity<WeightedNetwork> > L( Q );
		L.optimise();

		print_header(filename.c_str());
		L.print_basic();
		
	}
	


	return 0;
}
