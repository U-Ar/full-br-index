#include<iostream>
#include<cstdlib>

#include "br_index_fixed.hpp"

using namespace std;
using namespace bri;

long length = 0;

void help(){
	cout << "k_vs_space: breakdown of index space usage" << endl;
	cout << "Usage:       k_vs_space <length> <index>" << endl;
	cout << "   <index>   index file (with extension .bri/.bril)" << endl;
	exit(0);
}


int main(int argc, char** argv){

	if(argc < 2)
		help();

	br_index_fixed<> idx;
	cout << "Loading br-index length from" << argv[1] << endl;
	idx.load_from_file(argv[1]);
	cout << "--- Statistics of the text and the breakdown of the br-index space usage ---" << endl;

	for (ulint i = 0; i <= idx.border_length(); ++i)
	{
		idx.print_space(i);
	}

	cout << "------------------------------" << endl;
	cout << "full-directional br-index end." << endl;
	cout << "next: uni version of br-index with variable border length." << endl;
	cout << "------------------------------" << endl;

	for (ulint i = 0; i <= idx.border_length(); ++i)
	{
		idx.print_space_uni(i);
	}

}