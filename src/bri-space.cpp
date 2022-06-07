#include<iostream>

#include "br_index.hpp"

using namespace std;
using namespace bri;

bool nplcp = false;

void help(){
	cout << "bri-space: breakdown of index space usage" << endl;
	cout << "Usage:       bri-space <index>" << endl;
	cout << "   <index>   index file (with extension .bri)" << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc < 2)
		help();

	int ptr = 1;
	//while (ptr < argc-1) parse_args(argv, argc, ptr);

	br_index<> idx;

	cout << "Loading full-br-index" << endl;
	idx.load_from_file(argv[ptr]);
	cout << "--- Statistics of the text and the breakdown of the full-br-index space usage ---" << endl;
	
	auto space = idx.print_space();

}