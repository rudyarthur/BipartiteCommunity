#ifndef NETWORK_UTILS_H
#define NETWORK_UTILS_H

#include <fstream>

using namespace std; 

void print_header(const char* filename){

	ifstream infile(filename);
	if (infile.is_open() != true) {
		cerr << "The file " << filename << " does not exist" << endl;
		exit(1);
	}
  	
	string line;
	while(getline(infile, line)) {

		//find comments and blank lines: %, #, // allowed
		if(line.empty() || (line.find("%") == 0) || (line.find("#") == 0) || (line.find("//") == 0) ) { cout << line << endl;} 
			
	}
	infile.close();
	
}

#endif
