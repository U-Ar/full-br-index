#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index.hpp"
#include "br_index_nplcp.hpp"
#include "utils.hpp"

using namespace bri;
using namespace std;

string check = string();
long allowed = 0;
bool nplcp = false;
size_t left_len = 0;
size_t core_len = 0;

void help()
{
	cout << "bri-seedex: locate all occurrences of the input patterns" << endl;
    cout << "             with exact core and some mismatched characters."        << endl << endl;

	cout << "Usage: bri-seedex [options] <index> <patterns>" << endl;
    cout << "   -nplcp       use the version without PLCP." << endl;
    cout << "   -m <number>  max number of mismatched characters allowed (0 by default)" << endl;
	cout << "   -c <text>    check correctness of each pattern occurrence on this text file (must be the same indexed)" << endl;
	cout << "   <index>      index file (with extension .bri)" << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
    cout << "   <left>       length of the left region" << endl;
    cout << "   <core>       length of the core exact region" << endl;

	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if (s.compare("-c") == 0) 
    {

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -c option." << endl;
			help();
		}

		check = string(argv[ptr]);
		ptr++;
    
    }
    else if (s.compare("-m") == 0)
    {

        if(ptr>=argc-1){
            cout << "Error: missing parameter after -m option." << endl;
            help();
        }

        char* e;
        allowed = strtol(argv[ptr],&e,10);

        if(*e != '\0' || allowed < 0){
            cout << "Error: invalid negative value after -m option." << endl;
            help();
        }

        ptr++;

	}
    else if (s.compare("-nplcp") == 0)
    {

        nplcp = true;

    }
    else
    {

		cout << "Error: unknown option " << s << endl;
		help();

	}

}

template<class T>
void locate_all(ifstream& in, string patterns)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;

    string text;
    bool c = false;

    if (check.compare(string()) != 0)
    {
        c = true;

        ifstream ifs1(check);
        stringstream ss;
        ss << ifs1.rdbuf();
        text = ss.str();
    }

    auto t1 = high_resolution_clock::now();

    T idx;

    idx.load(in);

    auto t2 = high_resolution_clock::now();

    cout << "searching patterns with mismatches at most " << allowed << " ... " << endl;
    ifstream ifs(patterns);

    //read header of the pizza&chilli input file
	//header example:
	//# number=7 length=10 file=genome.fasta forbidden=\n\t
    string header;
    getline(ifs,header);

    ulint n = get_number_of_patterns(header);
    ulint m = get_patterns_length(header);

    ulint last_perc = 0;

    ulint occ_tot = 0;

    auto t3 = high_resolution_clock::now();
    auto t4 = high_resolution_clock::now();
    auto t5 = high_resolution_clock::now();
    ulint count_time = 0;
    ulint locate_time = 0;
    ulint tot_time = 0;

    // extract patterns from file and search them in the index
    for (ulint i = 0; i < n; ++i)
    {
        ulint perc = 100 * i / n;
        if (perc > last_perc)
        {
            cout << perc << "% done ..." << endl;
            last_perc = perc;
        }

        string p = string();

        for (ulint j = 0; j < m; ++j)
        {
            char c;
            ifs.get(c);
            p += c;
        }

        size_t m1 = left_len;
        size_t m2 = left_len + core_len;

        t3 = high_resolution_clock::now();
        auto samples = idx.seed_and_extend(p,m1,m2,allowed);
        t4 = high_resolution_clock::now();
        auto occs = idx.locate_samples(samples);
        t5 = high_resolution_clock::now();

        count_time += duration_cast<microseconds>(t4-t3).count();
        locate_time += duration_cast<microseconds>(t5-t4).count();
        occ_tot += occs.size();
        tot_time += duration_cast<microseconds>(t4-t3).count() + duration_cast<microseconds>(t5-t4).count();


        if (c) // check occurrences
        {
            cout << "number of occs with at most " << allowed << " mismatch   : " << occs.size() << endl;
            //cout << "the original pattern: " << p << endl;
            for (auto s : samples)
            {
                cout << s.second.range.first << " " << s.second.range.second << " " << s.second.j << " " << s.second.len << endl;
            }
            for (auto o : occs)
            {
                int mismatches = 0;
                for (size_t i = 0; i < m; ++i)
                {
                    if (text[o+i] != p[i]) mismatches++;
                }
                if (mismatches > allowed) 
                {
                    cout << "Error: wrong occurrence:  " << o << endl;
                    cout << "       original pattern:  " << p << endl;
                    cout << "       wrong    pattern:  " << text.substr(o,p.size()) << endl;
                }
                for (ulint k = m1; k < m2; ++k)
                {
                    if (text[o+k] != p[k]) 
                    {
                        cout << "Error: wrong occurrence:  " << o << endl;
                        cout << "       original pattern:  " << p << endl;
                        cout << "       wrong    pattern:  " << text.substr(o,p.size()) << endl;
                    }
                }
            }
        }

        
    }

    double occ_avg = (double)occ_tot / n;
    
    cout << endl << occ_avg << " average occurrences per pattern" << endl;

    ifs.close();

    ulint load = duration_cast<milliseconds>(t2-t1).count();
    cout << "Load time  : " << load << " milliseconds" << endl;

    cout << "Number of patterns n = " << n << endl;
	cout << "Pattern length     m = " << m << endl;
	cout << "Total number of occurrences   occt = " << occ_tot << endl << endl;

    cout << "LF-mapping time: " << count_time << " microseconds" << endl;
    cout << "Phi        time: " << locate_time << " microseconds" << endl;
    cout << "Total time : " << tot_time << " microseconds" << endl;
	cout << "Search time: " << (double)tot_time/n << " microseconds/pattern (total: " << n << " patterns)" << endl;
	cout << "Search time: " << (double)tot_time/occ_tot << " microseconds/occurrence (total: " << occ_tot << " occurrences)" << endl;
}



int main(int argc, char** argv)
{
    if (argc < 3) help();

    int ptr = 1;

    while (ptr < argc - 4) parse_args(argv, argc, ptr);

    string idx_file(argv[ptr]);
    string patt_file(argv[ptr+1]);

    char* e;
        
    left_len = strtol(argv[ptr+2],&e,10);

    if(*e != '\0'){
        cout << "Error: invalid value for <left>" << endl;
        help();
    }

    core_len = strtol(argv[ptr+3],&e,10);

    if(*e != '\0'){
        cout << "Error: invalid value for <core>" << endl;
        help();
    }

    ifstream in(idx_file);

    cout << "Loading br-index" << endl;

    if (nplcp)
        locate_all<br_index_nplcp<> >(in, patt_file);
    else 
        locate_all<br_index<> >(in, patt_file);

    in.close();

}