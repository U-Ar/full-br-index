#include<iostream>

#include "br_index.hpp"
#include "br_index_limited.hpp"

using namespace std;
using namespace bri;

bool limited = false;

void help(){
	cout << "bri-space: breakdown of index space usage" << endl;
	cout << "Usage:       bri-space [options] <index>" << endl;
	cout << "   -limited  use the limited version." << endl;
	cout << "   <index>   index file (with extension .bri/.bril)" << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if (s.compare("-limited") == 0)
    {

        limited = true;

    }
    else
    {
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv){

	if(argc < 2)
		help();

	int ptr = 1;
	while (ptr < argc-1) parse_args(argv, argc, ptr);

	if (!limited)
	{
		br_index<> idx;

		cout << "Loading full-br-index" << endl;
		idx.load_from_file(argv[ptr]);
		cout << "--- Statistics of the text and the breakdown of the br-index space usage ---" << endl;
		
		auto space = idx.print_space();
	}
	else 
	{
		br_index_limited<> idx;

		cout << "Loading br-index(limited)" << endl;
		idx.load_from_file(argv[ptr]);
		cout << "--- Statistics of the text and the breakdown of the br-index space usage ---" << endl;
		
		auto space = idx.print_space();
	}

}