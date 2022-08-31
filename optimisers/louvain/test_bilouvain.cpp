#include <iostream>
#include "../networks/network_utils.h"
#include "../optimiser.h"
#include "bilouvain.h"
#include <random>





using namespace std; 

int main(int argc, char **argv){

	vector< string > filenames = { "../../test_data/bipartite1.txt"};

	for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			
			BipartiteNetwork B(filename.c_str());		
			BarberModularity Q(B);
			print_header(filename.c_str());
			
			BiLouvain L( Q );
			L.optimise();
			L.print_basic();

	}

	return 0;

}
