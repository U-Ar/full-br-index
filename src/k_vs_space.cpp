#include<iostream>
#include<cstdlib>

#include "br_index_fixed.hpp"

using namespace std;
using namespace bri;

long length = 0;

void help(){
	cout << "k_vs_space: breakdown of index space usage" << endl;
	cout << "Usage:       k_vs_space <length> <index>" << endl;
	cout << "   <length>  use the limited version." << endl;
	cout << "   <index>   index file (with extension .bri/.bril)" << endl;
	exit(0);
}


int main(int argc, char** argv){

	if(argc < 3)
		help();

	char* e;
	length = strtol(argv[1],&e,10);

	if(*e != '\0' || length < 0){
		cout << "Error: invalid fixed length." << endl;
		return 0;
	}


	br_index_fixed<> idx;
	cout << "Loading br-index length " << length << " from" << argv[2] << endl;
	idx.load_from_file(argv[2]);
	cout << "--- Statistics of the text and the breakdown of the br-index space usage ---" << endl;

	for (ulint i = 0; i <= length; ++i)
	{
		idx.print_space(i);
	}

}