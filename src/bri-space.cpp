#include<iostream>

#include "br_index_full.hpp"

using namespace std;
using namespace bri;

bool limited = false;
bool fixed_length = false;

void help(){
	cout << "bri-space:   breakdown of index space usage" << endl;
	cout << "Usage:       bri-space [options] <index>" << endl;
	cout << "   <index>   index file (with extension .brif)" << endl;
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
	else if (s.compare("-fixed") == 0)
    {

        fixed_length = true;

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

	br_index_full<> idx;

	cout << "Loading fully functional br-index from " << argv[ptr] << " ..." << flush;
	idx.load_from_file(argv[ptr]);
	cout << "done."<< endl;
	cout << "--- The breakdown of the br-index space usage ---" << endl;
	auto space = idx.print_space();

	

}