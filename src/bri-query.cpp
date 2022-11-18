#include <unistd.h>

#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index.hpp"

using namespace bri;
using namespace std;

using hrclock=chrono::high_resolution_clock;
using chrono::duration_cast;
using chrono::duration;
using chrono::milliseconds;
using chrono::microseconds;

// struct containing command line parameters and other globals
struct Args {
    string check = "";
    int bl = -1;
    string query = "";
    string idx_file = "";
    string pattern_file = "";
    int k = -1;
    int t = -1;
};

void print_help(char** argv, Args &args) {
    std::cout << "Usage: " << argv[0] << " count  <index file> <pattern file> [options]" << std::endl;
    std::cout << "       " << argv[0] << " locate <index file> <pattern file> [options]" << std::endl;
    std::cout << "       " << argv[0] << " mems   <index file> <pattern file> [options]" << std::endl;
    std::cout << "       " << argv[0] << " full   <index file> <k(length)> <t(freq)> [options]" << std::endl;
    std::cout << "Compute searching query on given pattern using br-index." << std::endl;
    std::cout << "- count counts the number of the pattern occurring in the text." << std::endl;
    std::cout << "- locate locates all the starting positions of the pattern occurring in the text." << std::endl;
    std::cout << "- mems computes Maximal Exact Matches on the pattern and the text." << std::endl;
    std::cout << "- full computes full-task (total number of substrings of length at most k and occurring at least t times)" << std::endl;

    std::cout << std::endl << "  Options: " << std::endl
        << "\t-h  \tshow help and exit" << std::endl
        << "\t-l L\tparameter bl for contraction shortcut, def. index file's bl" << args.bl << std::endl;
    // << "\t-c C\tcheck correctness of each pattern occurrence" << std::endl
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
    if (argc == optind+4 && (strcmp(argv[optind],"full")==0)) {
        arg.query.assign(argv[optind]);
        arg.idx_file.assign(argv[optind+1]);
        sarg.assign( argv[optind+2] );
        arg.k = stoi(sarg);
        sarg.assign( argv[optind+3] );
        arg.t = stoi(sarg);
    }
    else if (argc == optind+3 && (strcmp(argv[optind],"count")==0 || strcmp(argv[optind],"locate")==0 || strcmp(argv[optind],"mems")==0 ))
    {
        arg.query.assign(argv[optind]);
        arg.idx_file.assign(argv[optind+1]);
        arg.pattern_file.assign(argv[optind+2]);
    } 
    else 
    {
        std::cout << "Invalid command line arguments" << std::endl;
        print_help(argv,arg);
    }
    // check algorithm parameters 
    if (arg.query=="full") {
        if (arg.k <= 0) {
            cout << "k must be positive integer" << endl;
            exit(1);
        }
        if (arg.t <= 0) {
            cout << "t must be positive integer" << endl;
            exit(1);
        }
    }
    //if(arg.bl < 0) {
    //    std::cout << "bl must be nonnegative integer\n";
    //    exit(1);
    //}
}

void query_count(Args& args, br_index& idx) {

    ifstream fp(args.pattern_file);
    if (!fp.is_open()) {
        cerr << "Cannot open pattern file: " << args.pattern_file << endl;
        exit(1);
    }
    string pattern((istreambuf_iterator<char>(fp)), istreambuf_iterator<char>());
    fp.close();

    cout << "Counting the pattern ... " << flush;

    auto t3 = hrclock::now();
    ulint cnt = idx.count(pattern);
    auto t4 = hrclock::now();

    cout << "done." << std::endl;

    cout << "#Occurrences: " << cnt << endl;
    cout << "Elapsed time: " << duration_cast<microseconds>(t4-t3).count() << "(us)" << endl;
}

void query_locate(Args& args, br_index& idx) {

    ifstream fp(args.pattern_file);
    if (!fp.is_open()) {
        cerr << "Cannot open pattern file: " << args.pattern_file << endl;
        exit(1);
    }
    string pattern((istreambuf_iterator<char>(fp)), istreambuf_iterator<char>());
    fp.close();

    cout << "Locating the pattern ... " << flush;

    auto t3 = hrclock::now();
    auto locations = idx.locate(pattern);
    auto t4 = hrclock::now();

    cout << "done." << std::endl;

    cout << "#Occurrences: " << locations.size() << endl;

    cout << "First min(#Occ,10) locations: ";
    ulint l = locations.size() >= 10 ? 10 : locations.size();
    for (size_t i = 0; i < l; ++i) cout << locations[i] << " ";
    cout << endl;

    cout << "Elapsed time: " << duration_cast<microseconds>(t4-t3).count() << "(us)" << endl;
}

void query_mems(Args& args, br_index& idx) {
    ifstream fp(args.pattern_file);
    if (!fp.is_open()) {
        cerr << "Cannot open pattern file: " << args.pattern_file << endl;
        exit(1);
    }
    string pattern((istreambuf_iterator<char>(fp)), istreambuf_iterator<char>());
    fp.close();

    cout << "MEMs on the pattern and the text ... " << flush;

    auto t3 = hrclock::now();
    ulint max_l = idx.maximal_exact_match(pattern);
    auto t4 = hrclock::now();

    cout << "done." << std::endl;

    cout << "Max length of maximal substrings: " << max_l << endl;

    cout << "Elapsed time: " << duration_cast<microseconds>(t4-t3).count() << "(us)" << endl;
}

void query_full(Args& args, br_index& idx) {
    ulint k = args.k, t = args.t;

    cout << "Max substring length k: " << k << endl;
    cout << "Min occurrence       t: " << t << endl; 
    cout << "Computing full-task ... " << flush;

    auto t3 = hrclock::now();
    ulint occ = idx.full_task(k,t);
    auto t4 = hrclock::now();

    cout << "done." << endl;

    cout << "#Occurrences: " << occ << endl;

    cout << "Elapsed time: " << duration_cast<microseconds>(t4-t3).count() << "(us)" << endl;
}

int main(int argc, char** argv)
{
    Args args;
    parse_args(argc, argv, args);

    auto t1 = hrclock::now();

    cout << "Loading br-index ... " << flush;

    br_index idx;
    ifstream fidx(args.idx_file);
    if (!fidx.is_open()) {
        cerr << "Cannot open index file: " << args.idx_file << endl;
        exit(1);
    }
    if (args.bl == -1) idx.load(fidx);
    else idx.load(fidx,args.bl);
    fidx.close();

    cout << "done." << endl;

    auto t2 = hrclock::now();

    cout << "Loading time: " << duration_cast<microseconds>(t2-t1).count() << "(us)" << endl;

    if (args.query == "full") {
        // string text((istreambuf_iterator<char>(fin)), istreambuf_iterator<char>());
        query_full(args,idx);
    }
    else if (args.query == "count") {
        query_count(args,idx);
    }
    else if (args.query == "locate") {
        query_locate(args,idx);
    } 
    else if (args.query == "mems") {
        query_mems(args,idx);
    }

}