#include <unistd.h>

#include <iostream>

#include "br_index.hpp"

using namespace std;
using namespace bri;

// struct containing command line parameters and other globals
struct Args {
    int bl = -1;
    string idx_file = "";
};

void print_help(char** argv, Args &args) {
    std::cout << "Usage: " << argv[0] << " <index file> [options]" << std::endl;
    std::cout << "Show breakdown of br-index space usage." << std::endl << std::endl;
    std::cout << "  Options: " << std::endl
        << "\t-h  \tshow help and exit" << std::endl 
        << "\t-l L\tparameter bl for contraction shortcut, def. index file's bl" << args.bl << std::endl;
    exit(1);
}

void parse_args( int argc, char** argv, Args& arg ) {
    int c;
    extern char *optarg;
    extern int optind;

    puts("==== Command line:");
    for(int i=0;i<argc;i++)
        printf(" %s",argv[i]);
    puts("");

    std::string sarg;
    while ((c = getopt( argc, argv, "l:h") ) != -1) {
        switch(c) {
            case 'l':
            sarg.assign( optarg );
            arg.bl = stoi( sarg ); break;
            case 'h':
            print_help(argv, arg); exit(1);
            case '?':
            cout << "Unknown option. Use -h for help." << endl;
            exit(1);
        }
    }
    if (argc == optind+1) {
        arg.idx_file.assign(argv[optind]);
    }
    else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_help(argv,arg);
    }
}

int main(int argc, char** argv) {

	Args arg;
    parse_args(argc, argv, arg);

	br_index idx;

	cout << "Loading br-index from " << arg.idx_file << " ..." << flush;

	ifstream in(arg.idx_file);
	if (arg.bl == -1) idx.load(in);
	else idx.load(in,arg.bl);

	cout << "done."<< endl;
	
	cout << "--- The breakdown of the br-index space usage ---" << endl;
	
	auto space = idx.print_space();


}