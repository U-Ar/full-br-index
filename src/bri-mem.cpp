#include <unistd.h>

#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index.hpp"

using namespace bri;
using namespace std;

// struct containing command line parameters and other globals
struct Args {
    string check = "";
    int bl = -1;
    string idx_file = "";
    string pattern_file = "";
};

void print_help(char** argv, Args &args) {
    std::cout << "Usage: " << argv[0] << " <index file> <pattern file> [options]" << std::endl;
    std::cout << "Compute Maximal Exact Matches on given pattern using br-index." << std::endl << std::endl;
    std::cout << "  Options: " << std::endl
        << "\t-c C\tcheck correctness of each pattern occurrence" << std::endl
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
    while ((c = getopt( argc, argv, "c:l:h") ) != -1) {
        switch(c) {
            case 'c':
            arg.check.assign(optarg); break;
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
    // the only input parameter is the file name 
    //if (argc == optind+1) {
    //    arg.input_file.assign( argv[optind] );
    //    if (!specify_out) arg.output_base.assign( argv[optind] );
    //}
    if (argc == optind+2) {
        arg.idx_file.assign(argv[optind]);
        arg.pattern_file.assign(argv[optind+1]);
    }
    else {
        std::cout << "Invalid number of arguments" << std::endl;
        print_help(argv,arg);
    }
    // check algorithm parameters 
    //if(arg.bl < 0) {
    //    std::cout << "bl must be nonnegative integer\n";
    //    exit(1);
    //}
}

template<class T>
void locate_all(Args& args, ifstream& in, string pattern_file)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;

    string text;
    bool c = false;

    if (args.check.compare(string()) != 0)
    {
        c = true;

        ifstream ifs1(args.check);
        stringstream ss;
        ss << ifs1.rdbuf();
        text = ss.str();
    }

    cout << "Loading br-index ... " << flush;

    auto t1 = high_resolution_clock::now();

    T idx;

    if (args.bl == -1) idx.load(in);
    else idx.load(in,args.bl);

    auto t2 = high_resolution_clock::now();

    cout << "done." << endl;

        
    ifstream ifs(args.pattern_file);

    stringstream ssp;
    ssp << ifs.rdbuf();
    string p = ssp.str();

    cout << "searching maximal exact matches ..." << endl;

    ulint m = p.size();

    ulint occ_tot = 0;

    auto t3 = high_resolution_clock::now();
    auto t4 = high_resolution_clock::now();
    ulint tot_time = 0;



    t3 = high_resolution_clock::now();

    idx.maximal_exact_match(p);

    t4 = high_resolution_clock::now();

    tot_time += duration_cast<microseconds>(t4-t3).count();


    /*if (c) // check occurrences
    {
        cout << "length of MEM : " << res.first << endl;
        cout << "number of MEM occs : " << occs.size() << endl;
        ulint maxlen = res.first;
        for (auto pair: res.second)
        {
            ulint offset = pair.first;
            auto occs = idx.locate_sample(pair.second);
            for (auto o: occs)
            {
                for (ulint i = 0; i < maxlen; ++i)
                {
                    if (text[o+i] != p[i+offset])
                    {
                        cout << "wrong occurrence  occ: " << o << " offset on P: " << offset << endl;
                    }
                }
            }
        }
    }*/

    

    ulint load = duration_cast<milliseconds>(t2-t1).count();
    cout << "Load time  : " << load << " milliseconds" << endl;

    cout << "Pattern length m = " << m << endl;
    
    // cout << "LF-mapping time: " << count_time << " microseconds" << endl;
    // cout << "Phi        time: " << locate_time << " microseconds" << endl;
    cout << "Computation time : " << tot_time << " microseconds" << endl;
    
}



int main(int argc, char** argv)
{
    Args arg;
    parse_args(argc, argv, arg);

    ifstream in(arg.idx_file);

    locate_all<br_index>(arg, in, arg.pattern_file);
    
    in.close();

}