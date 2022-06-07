#include<iostream>

#include "br_index.hpp"
#include "br_index_nplcp.hpp"

using namespace std;
using namespace bri;

bool nplcp = false;

void help(){
	cout << "bri-space: breakdown of index space usage" << endl;
	cout << "Usage:       bri-space [options] <index>" << endl;
	cout << "   -nplcp    use the version without PLCP." << endl;
	cout << "   <index>   index file (with extension .bri)" << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if (s.compare("-nplcp") == 0)
    {

        nplcp = true;

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

	if (nplcp)
	{
		br_index_nplcp<> idx;

		cout << "Loading br-index" << endl;
		idx.load_from_file(argv[ptr]);
		cout << "--- Statistics of the text and the breakdown of the br-index space usage ---" << endl;
		
		auto space = idx.print_space();
	}
	else 
	{
		br_index<> idx;

		cout << "Loading br-index" << endl;
		idx.load_from_file(argv[ptr]);
		cout << "--- Statistics of the text and the breakdown of the br-index space usage ---" << endl;
		
		auto space = idx.print_space();
	}

}